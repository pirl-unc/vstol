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


extern crate rayon;
extern crate serde;
use rayon::prelude::*;
use serde::{Serialize, Deserialize};
use std::collections::{HashMap, HashSet};
use crate::constants;
use crate::genomic_range::GenomicRange;
use crate::genomic_ranges_list::GenomicRangesList;
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

    /// Finds nearby variants.
    /// Please note that VariantCalls with the same super set of
    /// variant types will be identified as nearby variants.
    ///
    /// # Arguments
    /// * `query_variants_list`     -   VariantsList object.
    /// * `max_neighbor_distance`   -   maximum neighbor distance.
    /// * `num_threads`             -
    /// * `query_variant_types_map` -   HashMap where key is variant type and
    ///                                 value is a super set of the variant type.
    ///
    /// # Returns
    /// * `nearby_variants_map`     -   HashMap where key is variant ID and
    ///                                 value is a vector of variant IDs in
    ///                                 query_variants_list.
    pub fn find_nearby_variants(
        &self,
        query_variants_list: VariantsList,
        max_neighbor_distance: isize,
        num_threads: usize,
        query_variant_types_map: &HashMap<&str, String>) -> HashMap<String, Vec<String>> {
        // Step 1. Split Variant objects by (chromosome_1, chromosome_2, variant_type)
        let mut variants_map: HashMap<(String, String, String), Vec<Variant>> = HashMap::new();
        for variant in self.variants.iter() {
            // Get the query variant type
            let query_variant_type: String = query_variant_types_map.get(variant.variant_calls[0].variant_type.as_str()).unwrap().to_string();
            let key = (
                variant.variant_calls[0].chromosome_1.clone(),
                variant.variant_calls[0].chromosome_2.clone(),
                query_variant_type.clone(),
            );
            variants_map
                .entry(key)
                .or_insert(Vec::new())
                .push(variant.clone());
        }

        // Step 2. Split query Variant objects by (chromosome_1, chromosome_2, variant_type)
        let mut query_variants_map: HashMap<(String, String, String), Vec<Variant>> = HashMap::new();
        for variant in query_variants_list.variants.iter() {
            // Get the query variant type
            let query_variant_type: String = query_variant_types_map.get(variant.variant_calls[0].variant_type.as_str()).unwrap().to_string();
            let key = (
                variant.variant_calls[0].chromosome_1.clone(),
                variant.variant_calls[0].chromosome_2.clone(),
                query_variant_type.clone(),
            );
            query_variants_map
                .entry(key)
                .or_insert(Vec::new())
                .push(variant.clone());
        }

        // Step 3. Find nearby variants
        let thread_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();
        let results: Vec<HashMap<String, Vec<String>>> = thread_pool.install(|| {
            variants_map.par_iter().map(|(key, values)| {
                // Clone Variant objects
                let mut variants_: Vec<Variant> = values.clone();
                let mut nearby_variants_map_: HashMap<String, Vec<String>> = HashMap::new(); // key = Variant.id, value = Vec<Query Variant.id>
                if query_variants_map.contains_key(key) {
                    let mut query_variants_: Vec<Variant> = query_variants_map.get(key).unwrap().clone();
                    let mut max_neighbor_distance_: isize = max_neighbor_distance;
                    if variants_[0].variant_calls[0].variant_type == constants::SINGLE_NUCLEOTIDE_VARIANT ||
                        variants_[0].variant_calls[0].variant_type == constants::MULTI_NUCLEOTIDE_VARIANT {
                        max_neighbor_distance_ = 0;
                    }

                    // Sort Variant objects
                    variants_.sort_by(|a, b| a.variant_calls[0].position_1.cmp(&b.variant_calls[0].position_1));
                    query_variants_.sort_by(|a, b| a.variant_calls[0].position_1.cmp(&b.variant_calls[0].position_1));

                    // Find nearby query Variant objects
                    for variant_ in &variants_ {
                        for query_variant_ in &query_variants_ {
                            let mut include: bool = false;
                            for variant_call_ in &variant_.variant_calls {
                                for query_variant_call_ in &query_variant_.variant_calls {
                                    let distance_1: isize = (variant_call_.position_1 - query_variant_call_.position_1).abs();
                                    let distance_2: isize = (variant_call_.position_2 - query_variant_call_.position_2).abs();
                                    if (distance_1 <= max_neighbor_distance_) && (distance_2 <= max_neighbor_distance_) {
                                        include = true;
                                    }
                                }
                            }
                            if include {
                                nearby_variants_map_
                                    .entry(variant_.id.to_string())
                                    .or_insert(Vec::new())
                                    .push(query_variant_.id.to_string());
                            }
                        }
                    }
                }
                nearby_variants_map_
            })
            .collect()
        });

        // Step 4. Collect IDs into one HashMap
        let mut nearby_variants_map: HashMap<String, Vec<String>> = HashMap::new(); // key = Variant.id, value = Vec<Query Variant.id>
        for result in &results {
            for (key, value) in result.iter() {
                for variant_id in value {
                    nearby_variants_map
                        .entry(key.to_string())
                        .or_insert(Vec::new())
                        .push(variant_id.to_string());
                }
            }
        }

        return nearby_variants_map;
    }

    /// Finds variants overlapping query regions.
    ///
    /// # Arguments
    /// * `genomic_regions_list`    -   VariantsList object.
    /// * `padding`                 -   padding to apply to GenomicRange start and end.
    /// * `num_threads`             -   number of threads.
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

    /// This function returns a VariantsList object comprised of Variant objects
    /// that have intersecting VariantCall objects across all supplied
    /// VariantsList objects. Note that VariantCalls with the same super set of
    /// variant types will be merged into one.
    ///
    /// # Arguments
    /// * `variants_lists`              -   vector of VariantsList objects.
    /// * `num_threads`                 -   number of threads.
    /// * `max_neighbor_distance`       -   maximum neighbor distance.
    /// * `query_variant_types_map`     -   HashMap where key is variant type and
    ///                                     value is a super set of the variant type.
    ///
    /// # Returns
    /// * `intersecting_variants_list`  -   intersecting VariantsList object.
    pub fn intersect(
        variants_lists: Vec<VariantsList>,
        num_threads: usize,
        max_neighbor_distance: isize,
        query_variant_types_map: &HashMap<&str, String>) -> VariantsList {
        // Step 1. Split variants by (chromosome_1, chromosome_2, super variant_type)
        let mut variant_calls_map: HashMap<(String, String, String), Vec<VariantCall>> = HashMap::new();
        let mut variants_list_id: isize = 0;
        for variants_list in variants_lists.iter() {
            for variant in variants_list.variants.iter() {
                for variant_call in variant.variant_calls.iter() {
                    let mut variant_call_ = variant_call.clone();

                    // Denote which variants list this VariantCall is from
                    variant_call_.add_attribute("variants_list_id".to_string(), variants_list_id.to_string());

                    // Get the query variant type
                    let query_variant_type: String = query_variant_types_map.get(variant_call_.variant_type.as_str()).unwrap().to_string();
                    let key = (
                        variant_call_.chromosome_1.clone(),
                        variant_call_.chromosome_2.clone(),
                        query_variant_type.clone(),
                    );
                    variant_calls_map
                        .entry(key)
                        .or_insert(Vec::new())
                        .push(variant_call_);
                }
            }
            variants_list_id += 1;
        }

        // Step 2. Assign variant IDs
        let thread_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();
        let results: Vec<HashMap<String, String>> = thread_pool.install(|| {
            variant_calls_map.par_iter().map(|(key, values)| {
                // Clone values
                let mut variant_calls_temp: Vec<VariantCall> = values.to_vec();

                // Sort variant calls
                variant_calls_temp.sort_by(|a, b| a.position_1.cmp(&b.position_1));

                // Maximum neighbor distance must be 0 for SNVs and MNVs
                let mut max_neighbor_distance_: isize = max_neighbor_distance;
                if variant_calls_temp[0].variant_type == constants::SINGLE_NUCLEOTIDE_VARIANT ||
                    variant_calls_temp[0].variant_type == constants::MULTI_NUCLEOTIDE_VARIANT {
                    max_neighbor_distance_ = 0;
                }

                // Iterate through the vector of VariantCall objects,
                // which all have the same
                // (chromosome_1, chromosome_2, super variant_type),
                // and assign variant IDs
                let mut variant_idx = 1;
                let mut variant_ids_map: HashMap<String, String> = HashMap::new(); // key = VariantCall.id, value = Variant.id
                for i in 0..variant_calls_temp.len() {
                    if !variant_ids_map.contains_key(&variant_calls_temp[i].id) {
                        // Iterate through the subsequent (sorted) VariantCall objects
                        // and figure out the same variant calls as long as the
                        // neighbor distance is within the desired maximum distance
                        let mut matched_variant_call_ids: Vec<String> = Vec::new();
                        for j in (i + 1)..variant_calls_temp.len() {
                            let distance_1: isize = (variant_calls_temp[i].position_1 - variant_calls_temp[j].position_1).abs();
                            let distance_2: isize = (variant_calls_temp[i].position_2 - variant_calls_temp[j].position_2).abs();
                            if (distance_1 <= max_neighbor_distance_) && (distance_2 <= max_neighbor_distance_) {
                                if let (Some(vl_id_1), Some(vl_id_2)) = (variant_calls_temp[i].attributes.get("variants_list_id"), variant_calls_temp[j].attributes.get("variants_list_id")) {
                                    if vl_id_1 != vl_id_2 {
                                        // These VariantCall objects are from different VariantsList objects
                                        matched_variant_call_ids.push(variant_calls_temp[j].id.clone());
                                    }
                                }
                            } else {
                                if distance_1 > max_neighbor_distance_ {
                                    break;
                                }
                            }
                        }
                        if !matched_variant_call_ids.is_empty() {
                            // <chromosome_1>-<chromosome_2>-<variant_type>-<variant_id>
                            let variant_id = format!("{}-{}-{}_{}", key.0, key.1, key.2, variant_idx.to_string());
                            // Add the current VariantCall ID into variant_ids_map
                            variant_ids_map.insert(variant_calls_temp[i].id.clone(), variant_id.clone());
                            // Add all the matched VariantCall IDs into variant_ids_map
                            for j in 0..matched_variant_call_ids.len() {
                                variant_ids_map.insert(matched_variant_call_ids[j].clone(), variant_id.clone());
                            }
                            variant_idx += 1;
                        }
                    }
                }
                variant_ids_map
            })
            .collect()
        });

        // Step 3. Collect all VariantCall and Variant IDs into one HashMap
        let mut variant_ids_map: HashMap<String, String> = HashMap::new(); // key = VariantCall.id, value = Variant.id
        for result in &results {
            for (key, value) in result.iter() {
                variant_ids_map.insert(key.to_string(), value.to_string());
            }
        }

        // Step 4. Merge VariantCall objects into common Variant objects
        let mut variants_map: HashMap<String, Variant> = HashMap::new();
        for variants_list in variants_lists.iter() {
            for variant in variants_list.variants.iter() {
                for variant_call in variant.variant_calls.iter() {
                    if !variant_ids_map.contains_key(&variant_call.id) {
                        // Current VariantCall does not intersect with any
                        // other VariantCall objects
                        continue;
                    }
                    // Retrieve variant ID
                    let variant_id = variant_ids_map
                        .get(&variant_call.id)
                        .map(|s| s.to_string())
                        .unwrap_or_else(|| String::new());
                    if let Some(matched_variant) = variants_map.get_mut(&variant_id) {
                        matched_variant.add_variant_call(variant_call.clone());
                    } else {
                        let mut new_variant = Variant::new(variant_id.clone());
                        new_variant.add_variant_call(variant_call.clone());
                        variants_map.insert(variant_id.clone(), new_variant);
                    }
                }
            }
        }

        // Step 5. Create a merged VariantsList object
        let mut variants_list_merged = VariantsList::new();
        for (_variant_id, variant) in variants_map.iter() {
            variants_list_merged.add_variant(variant.clone());
        }

        return variants_list_merged;
    }

    /// This function merges a vector of VariantsList objects into one.
    /// Please note that VariantCalls with the same super set of
    /// variant types will be merged into one.
    ///
    /// # Arguments
    /// * `variants_lists`          -   vector of VariantsList objects.
    /// * `num_threads`             -   number of threads.
    /// * `max_neighbor_distance`   -   maximum neighbor distance.
    /// * `query_variant_types_map` -   HashMap where key is variant type and
    ///                                 value is a super set of the variant type.
    ///
    /// # Returns
    /// * `merged_variants_list`    -   merged VariantsList object.
    pub fn merge(
        variants_lists: Vec<VariantsList>,
        num_threads: usize,
        max_neighbor_distance: isize,
        query_variant_types_map: &HashMap<&str, String>) -> VariantsList {
        // Step 1. Split variants by (chromosome_1, chromosome_2, super variant_type)
        let mut variant_calls_map: HashMap<(String, String, String), Vec<VariantCall>> = HashMap::new();
        let mut variants_list_id: isize = 0;
        for variants_list in variants_lists.iter() {
            for variant in variants_list.variants.iter() {
                for variant_call in variant.variant_calls.iter() {
                    let mut variant_call_ = variant_call.clone();

                    // Denote which variants list this VariantCall is from
                    variant_call_.add_attribute("variants_list_id".to_string(), variants_list_id.to_string());

                    // Get the query variant type
                    let query_variant_type: String = query_variant_types_map.get(variant_call_.variant_type.as_str()).unwrap().to_string();
                    let key = (
                        variant_call_.chromosome_1.clone(),
                        variant_call_.chromosome_2.clone(),
                        query_variant_type.clone()
                    );
                    variant_calls_map
                        .entry(key)
                        .or_insert(Vec::new())
                        .push(variant_call_);
                }
            }
            variants_list_id += 1;
        }

        // Step 2. Assign variant IDs
        let thread_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();
        let results: Vec<HashMap<String, String>> = thread_pool.install(|| {
            variant_calls_map.par_iter().map(|(key, values)| {
                // Clone values
                let mut variant_calls_temp: Vec<VariantCall> = values.to_vec();

                // Sort variant calls
                variant_calls_temp.sort_by(|a, b| a.position_1.cmp(&b.position_1));

                // Maximum neighbor distance must be 0 for SNVs and MNVs
                let mut max_neighbor_distance_: isize = max_neighbor_distance;
                if variant_calls_temp[0].variant_type == constants::SINGLE_NUCLEOTIDE_VARIANT ||
                    variant_calls_temp[0].variant_type == constants::MULTI_NUCLEOTIDE_VARIANT {
                    max_neighbor_distance_ = 0;
                }

                // Iterate through the vector of VariantCall objects,
                // which all have the same
                // (chromosome_1, chromosome_2, super variant_type),
                // and assign variant IDs
                let mut variant_idx = 1;
                let mut variant_ids_map: HashMap<String, String> = HashMap::new(); // key = VariantCall.id, value = Variant.id
                for i in 0..variant_calls_temp.len() {
                    if !variant_ids_map.contains_key(&variant_calls_temp[i].id) {
                        // <chromosome_1>-<chromosome_2>-<variant_type>-<variant_id>
                        let variant_id = format!("{}-{}-{}_{}", key.0, key.1, key.2, variant_idx.to_string());
                        variant_ids_map.insert(variant_calls_temp[i].id.clone(), variant_id.clone());
                        // Iterate through the subsequent (sorted) VariantCall objects
                        // and assign the same variant ID as long as the
                        // neighbor distance is within the desired maximum distance
                        for j in (i + 1)..variant_calls_temp.len() {
                            let distance_1: isize = (variant_calls_temp[i].position_1 - variant_calls_temp[j].position_1).abs();
                            let distance_2: isize = (variant_calls_temp[i].position_2 - variant_calls_temp[j].position_2).abs();
                            if (distance_1 <= max_neighbor_distance_) && (distance_2 <= max_neighbor_distance_) {
                                variant_ids_map.insert(variant_calls_temp[j].id.clone(), variant_id.clone());
                            } else {
                                if distance_1 > max_neighbor_distance_ {
                                    break;
                                }
                            }
                        }
                        variant_idx += 1;
                    }
                }
                variant_ids_map
            })
            .collect()
        });

        // Step 3. Collect all VariantCall and Variant IDs into one HashMap
        let mut variant_ids_map: HashMap<String, String> = HashMap::new(); // key = VariantCall.id, value = Variant.id
        for result in &results {
            for (key, value) in result.iter() {
                variant_ids_map.insert(key.to_string(), value.to_string());
            }
        }

        // Step 4. Merge VariantCall objects into common Variant objects
        let mut variants_map: HashMap<String, Variant> = HashMap::new();
        for variants_list in variants_lists.iter() {
            for variant in variants_list.variants.iter() {
                for variant_call in variant.variant_calls.iter() {
                    let variant_id = variant_ids_map
                        .get(&variant_call.id)
                        .map(|s| s.to_string())
                        .unwrap_or_else(|| String::new());
                    if let Some(matched_variant) = variants_map.get_mut(&variant_id) {
                        matched_variant.add_variant_call(variant_call.clone());
                    } else {
                        let mut new_variant = Variant::new(variant_id.clone());
                        new_variant.add_variant_call(variant_call.clone());
                        variants_map.insert(variant_id.clone(), new_variant);
                    }
                }
            }
        }

        // Step 5. Create a merged VariantsList object
        let mut variants_list_merged = VariantsList::new();
        for (_variant_id, variant) in variants_map.iter() {
            variants_list_merged.add_variant(variant.clone());
        }

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

