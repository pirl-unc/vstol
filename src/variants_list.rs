// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


extern crate log;
extern crate rayon;
extern crate serde;
use log::info;
use rayon::prelude::*;
use serde::{Serialize, Deserialize};
use std::collections::{HashMap, HashSet};
use std::process;
use crate::constants;
use crate::genomic_range::GenomicRange;
use crate::genomic_ranges_list::GenomicRangesList;
use crate::utilities::find_clusters;
use crate::variant::Variant;
use crate::variant_call::VariantCall;
use crate::variant_filter::VariantFilter;


#[derive(Debug, Serialize, Deserialize)]
pub struct VariantsList {
    pub variants: Vec<Variant>
}

impl VariantsList {
    pub fn new() -> Self {
        Self {
            variants: Vec::new()
        }
    }

    pub fn add_variant(&mut self, variant: Variant) {
        self.variants.push(variant);
    }

    pub fn filter(&self, variant_filters: Vec<VariantFilter>, num_threads: usize) -> VariantsList {
        // Step 1. Identify variant IDs to keep
        let thread_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();
        let variant_ids_to_keep: Vec<String> = thread_pool.install(|| {
            self.variants.par_iter().map(|variant| {
                let mut keep = true;
                for variant_filter in &variant_filters {
                    if variant_filter.keep(&variant) == false {
                        keep = false;
                        break;
                    }
                }
                if keep {
                    return variant.id.to_string();
                } else {
                    return "".to_string();
                }
            }).collect()
        });

        // Step 2. Create a HashSet of variant IDs to keep
        let variant_ids_to_keep_set: HashSet<String> = variant_ids_to_keep.into_iter().collect();

        // Step 3. Create a filtered VariantsList
        let mut filtered_variants_list: VariantsList = VariantsList::new();
        for variant in &self.variants {
            if variant_ids_to_keep_set.contains(&variant.id) {
                filtered_variants_list.add_variant(variant.clone());
            }
        }

        return filtered_variants_list;
    }

    /// Finds variants overlapping query regions.
    ///
    /// # Arguments
    /// * `genomic_regions_list`    -   VariantsList object.
    /// * `padding`                 -   Padding to apply to GenomicRange start and end.
    /// * `num_threads`             -   Number of threads.
    ///
    /// # Returns
    /// * `nearby_variants_map`     -   HashMap where key is variant call ID and
    ///                                 value is a vector of GenomicRange IDs.
    pub fn overlap(
        &self,
        genomic_regions_list: GenomicRangesList,
        num_threads: usize,
        padding: isize) -> HashMap<String, Vec<String>> {
         // Step 1. Split Variant objects by chromosome
        let mut variants_map: HashMap<String, Vec<Variant>> = HashMap::new();
        for variant in self.variants.iter() {
            let key1 = variant.variant_calls[0].chromosome_1.clone();
            variants_map
                .entry(key1)
                .or_insert(Vec::new())
                .push(variant.clone());
            if variant.variant_calls[0].chromosome_1 != variant.variant_calls[0].chromosome_2 {
                let key2 = variant.variant_calls[0].chromosome_2.clone();
                variants_map
                    .entry(key2)
                    .or_insert(Vec::new())
                    .push(variant.clone());
            }
        }

        // Step 2. Filter regions
        let thread_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();
        let results: Vec<HashMap<String, Vec<String>>> = thread_pool.install(|| {
            variants_map.par_iter().map(|(key, variants)| {
                let mut result: HashMap<String, Vec<String>> = HashMap::new();
                if genomic_regions_list.genomic_ranges_map.contains_key(key) {
                    let genomic_ranges: Vec<GenomicRange> = genomic_regions_list.genomic_ranges_map.get(key).unwrap().clone();
                    for variant in variants {
                        for variant_call in &variant.variant_calls {
                            for genomic_range in &genomic_ranges {
                                let start_1 = variant_call.position_1 - padding;
                                let end_1 = variant_call.position_1 + padding;
                                let start_2 = variant_call.position_2 - padding;
                                let end_2 = variant_call.position_2 + padding;
                                let overlaps_1 = genomic_range.overlaps(variant_call.chromosome_1.to_string(), start_1, end_1);
                                let overlaps_2 = genomic_range.overlaps(variant_call.chromosome_2.to_string(), start_2, end_2);
                                if overlaps_1 || overlaps_2 {
                                    result
                                        .entry(variant_call.id.clone())
                                        .or_insert(Vec::new())
                                        .push(genomic_range.id().clone());
                                }
                            }
                        }
                    }
                }
                result
            })
            .collect()
        });

        // Step 3. Collect all results into one HashMap
        let mut overlap_regions_map: HashMap<String, Vec<String>> = HashMap::new();
        for result in &results {
            for (key, values) in result.iter() {
                for value in values {
                    overlap_regions_map
                        .entry(key.clone())
                        .or_insert(Vec::new())
                        .push(value.clone());
                }
            }
        }

        return overlap_regions_map;
    }

    /// This function merges a vector of VariantsList objects into one.
    ///
    /// # Arguments
    /// * `variants_lists`              -   Vector of VariantsList objects.
    /// * `num_threads`                 -   Number of threads.
    /// * `max_neighbor_distance`       -   Maximum neighbor distance.
    /// * `match_all_breakpoints`       -   If true, both pairs of breakpoints of two variant calls
    ///                                     must be near each other (start1 == start2 && end1 == end2).
    ///                                     If false, only one pair of breakpoints of two variant calls
    ///                                     must be near each other (start1 == start2 || end1 == end2).
    /// * `match_variant_types`         -   If true, variant types (super types) must match.
    /// * `intersect_only`              -   If true, private variant calls are excluded (i.e. should be true for intersect).
    /// * `variant_types_map`           -   HashMap where key is variant type and
    ///                                     value is a super set of the variant type.
    ///
    /// # Returns
    /// * `merged_variants_list`        -   Merged VariantsList object.
    pub fn merge(
        variants_lists: Vec<VariantsList>,
        num_threads: usize,
        max_neighbor_distance: isize,
        match_all_breakpoints: bool,
        match_variant_types: bool,
        intersect_only: bool,
        variant_types_map: &HashMap<&str, String>) -> VariantsList {
        let thread_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();

        // Step 1. Split variant calls by breakpoint chromosome
        // key      = chromosome
        // value    = Vec<(position,VariantCall)>
        info!("Started splitting VariantCall objects by breakpoint chromosome");
        let mut variant_calls_map: HashMap<String, Vec<(isize,VariantCall)>> = HashMap::new();
        let mut variants_list_id: isize = 0;
        for variants_list in variants_lists.iter() {
            for variant in variants_list.variants.iter() {
                for variant_call in variant.variant_calls.iter() {
                    let mut variant_call_ = variant_call.clone();

                    // Save VariantsList ID and Variant ID
                    variant_call_.add_attribute(
                        "temp_variants_list_id".to_string(),
                        variants_list_id.to_string()
                    );
                    variant_call_.add_attribute(
                        "temp_variant_id".to_string(),
                        variant.id.to_string()
                    );

                    // Append chromosome_1
                    variant_calls_map
                        .entry(variant_call_.chromosome_1.to_string())
                        .or_insert(Vec::new())
                        .push((variant_call_.position_1, variant_call_.clone()));

                    // Append chromosome_2
                    variant_calls_map
                        .entry(variant_call_.chromosome_2.to_string())
                        .or_insert(Vec::new())
                        .push((variant_call_.position_2, variant_call_.clone()));
                }
            }
            variants_list_id += 1;
        }
        info!("Finished splitting VariantCall objects by breakpoint chromosome");

        // Step 2. Sort variant_records_map by position
        info!("Started sorting VariantCall objects");
        thread_pool.install(|| {
            variant_calls_map.par_iter_mut().for_each(|(_key, variant_call_positions)| {
                variant_call_positions.sort_by_key(|&(position,_)| position);
            });
        });
        info!("Finished sorting VariantCall objects");

        // Step 3. Identify nearby variant calls
        info!("Started identifying nearby VariantCall objects");
        let results: HashSet<(String,String)> = thread_pool.install(|| {
            variant_calls_map.par_iter().flat_map(|(chromosome, variant_call_positions)| {
                // Check if breakpoints intersect
                let mut pairs: HashSet<(String,String)> = HashSet::new(); // (VariantCall.id, VariantCall.id)
                for i in 0..variant_call_positions.len() {
                    let position_1: isize = variant_call_positions[i].0;
                    let variant_call_1: &VariantCall = &variant_call_positions[i].1;
                    let variant_type_1: &str = variant_types_map.get(variant_call_1.variant_type.as_str()).unwrap().as_str();
                    for j in (i + 1)..variant_call_positions.len() {
                        let position_2: isize = variant_call_positions[j].0;
                        let variant_call_2: &VariantCall = &variant_call_positions[j].1;
                        let variant_type_2: &str = variant_types_map.get(variant_call_2.variant_type.as_str()).unwrap().as_str();
                        let distance: isize = (position_1 - position_2).abs();
                        if (distance > max_neighbor_distance) {
                            break;
                        }
                        if variant_call_1.id == variant_call_2.id {
                            continue;
                        }
                        let mut cluster: bool = false;
                        if let (Some(vl_id_1), Some(vl_id_2), Some(v_id_1), Some(v_id_2)) = (
                            variant_call_1.attributes.get("temp_variants_list_id"),
                            variant_call_2.attributes.get("temp_variants_list_id"),
                            variant_call_1.attributes.get("temp_variant_id"),
                            variant_call_2.attributes.get("temp_variant_id")
                        ) {
                            // The VariantCall objects are from the same list
                            // and were originally in the same Variant
                            if (vl_id_1 == vl_id_2) && (v_id_1 == v_id_2) && (variant_call_1.id != variant_call_2.id) {
                                cluster = true;
                            }

                            // The VariantCall objects are from different lists
                            if (vl_id_1 != vl_id_2) {
                                if (match_variant_types == false) || ((match_variant_types == true) && (variant_type_1 == variant_type_2)) {
                                    // Maximum neighbor distance must be 0 for SNVs and MNVs
                                    let mut max_neighbor_distance_: isize = max_neighbor_distance;
                                    if (variant_call_1.variant_type == constants::SINGLE_NUCLEOTIDE_VARIANT) ||
                                        (variant_call_1.variant_type == constants::MULTI_NUCLEOTIDE_VARIANT) {
                                        // variant_call_1 and variant_call_2 are both SNV or MNV
                                        max_neighbor_distance_ = 0;
                                    }

                                    let distance_11: isize = (variant_call_1.position_1 - variant_call_2.position_1).abs();
                                    let distance_22: isize = (variant_call_1.position_2 - variant_call_2.position_2).abs();
                                    let distance_12: isize = (variant_call_1.position_1 - variant_call_2.position_2).abs();
                                    let distance_21: isize = (variant_call_1.position_2 - variant_call_2.position_1).abs();
                                    if (match_all_breakpoints == true) {
                                        if ((variant_call_1.chromosome_1 == variant_call_2.chromosome_1) &&
                                            (variant_call_1.chromosome_2 == variant_call_2.chromosome_2) &&
                                            (distance_11 <= max_neighbor_distance_) &&
                                            (distance_22 <= max_neighbor_distance_)) ||
                                           ((variant_call_1.chromosome_1 == variant_call_2.chromosome_2) &&
                                            (variant_call_1.chromosome_2 == variant_call_2.chromosome_1) &&
                                            (distance_12 <= max_neighbor_distance_) &&
                                            (distance_21 <= max_neighbor_distance_)) {
                                            cluster = true;
                                        }
                                    } else {
                                        if ((variant_call_1.chromosome_1 == variant_call_2.chromosome_1) && (distance_11 <= max_neighbor_distance_)) ||
                                            ((variant_call_1.chromosome_2 == variant_call_2.chromosome_2) && (distance_22 <= max_neighbor_distance_)) ||
                                            ((variant_call_1.chromosome_1 == variant_call_2.chromosome_2) && (distance_12 <= max_neighbor_distance_)) ||
                                            ((variant_call_1.chromosome_2 == variant_call_2.chromosome_1) && (distance_21 <= max_neighbor_distance_)) {
                                            cluster = true;
                                        }
                                    }
                                }
                            }

                            if cluster {
                                if (pairs.contains(&(variant_call_1.id.to_string(), variant_call_2.id.to_string())) == false) &&
                                    (pairs.contains(&(variant_call_2.id.to_string(), variant_call_1.id.to_string())) == false) {
                                    pairs.insert((variant_call_1.id.to_string(), variant_call_2.id.to_string()));
                                }
                            }
                        }
                    }
                }
                pairs
            })
            .collect()
        });
        info!("Finished identifying nearby VariantCall objects");

        // Step 4. By transitive property, identify clusters
        info!("Started identifying clusters");
        let clusters = find_clusters(results);
        info!("Finished identifying clusters");

        // Step 5. Assign Variant ID to each VariantCall
        // key      =   VariantCall.id
        // value    =   Varinat.id
        info!("Started mapping VariantCall to Variant");
        let mut variant_ids_map: HashMap<String, String> = HashMap::new();
        let mut new_variant_id = 1;
        for variant_call_ids in clusters.iter() {
            for variant_call_id in variant_call_ids.iter() {
                variant_ids_map.insert(variant_call_id.to_string(), new_variant_id.to_string());
            }
            new_variant_id += 1;
        }
        info!("Finished mapping VariantCall to Variant");

        // Step 6. Merge VariantCall objects into common Variant objects
        info!("Started merging VariantCall objects into common Variant objects");
        let mut variants_map: HashMap<String, Variant> = HashMap::new();
        for variants_list in variants_lists.iter() {
            for variant in variants_list.variants.iter() {
                for variant_call in variant.variant_calls.iter() {
                    if variant_ids_map.contains_key(&variant_call.id) {
                        // Retrieve variant ID
                        let variant_id: String = variant_ids_map
                            .get(&variant_call.id)
                            .map(|s| s.to_string())
                            .unwrap_or_else(|| String::new());

                        // Fetch the corresponding Variant
                        if variants_map.contains_key(&variant_id) {
                            variants_map.entry(variant_id.clone()).and_modify(|v| {
                                v.add_variant_call(variant_call.clone());
                            });
                        } else {
                            let mut new_variant = Variant::new(variant_id.clone());
                            new_variant.add_variant_call(variant_call.clone());
                            variants_map.insert(new_variant.id.clone(), new_variant);
                        }
                    } else {
                        if intersect_only == false {
                            // Current VariantCall does not intersect with any
                            // other VariantCall objects
                            let mut new_variant = Variant::new(new_variant_id.to_string());
                            new_variant.add_variant_call(variant_call.clone());
                            variants_map.insert(new_variant.id.clone(), new_variant);
                            new_variant_id += 1;
                        }
                    }
                }
            }
        }
        info!("Finished merging VariantCall objects into common Variant objects");

        // Step 7. Create a merged VariantsList object
        info!("Started merging Variant objects into a common merged VariantsList");
        let mut variants_list_merged = VariantsList::new();
        for (_variant_id, variant) in variants_map.iter() {
            variants_list_merged.add_variant(variant.clone());
        }
        info!("Finished merging Variant objects into a common merged VariantsList");
        return variants_list_merged;
    }
}

impl Clone for VariantsList {
    fn clone(&self) -> Self {
        let mut variants: Vec<Variant> = Vec::new();
        for variant in &self.variants {
            variants.push(variant.clone());
        }
        VariantsList {
            variants: variants
        }
    }
}

