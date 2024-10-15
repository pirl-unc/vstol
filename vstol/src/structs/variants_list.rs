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
use rayon::prelude::*;
use serde::{Serialize, Deserialize};
use std::cmp::{max,min};
use std::collections::{HashMap, HashSet};
use crate::constants;
use crate::algorithms::clustering::find_clusters;
use crate::structs::genomic_range::GenomicRange;
use crate::structs::genomic_ranges_list::GenomicRangesList;
use crate::structs::variant::Variant;
use crate::structs::variant_call::VariantCall;
use crate::structs::variant_filter::VariantFilter;


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

    /// Add a Variant object.
    pub fn add_variant(&mut self, variant: Variant) {
        self.variants.push(variant);
    }

    /// Compare two VariantsList objects and return a tuple of three VariantsList:
    /// 1. Shared between A and B.
    /// 2. Specific to A.
    /// 3. Specific to B.
    ///
    /// # Arguments
    ///
    /// * If `match_all_breakpoints==true`, both pairs of breakpoints of two variant calls must be
    /// near each other `(start1 == start2 && end1 == end2)`.
    ///
    /// * If `match_all_breakpoints==false`, only one pair of breakpoints of two variant calls
    /// must be near each other `(start1 == start2 || end1 == end2)`.
    ///
    /// * If `match_variant_types==true`, variant types (super types) must match.
    ///
    /// * `min_ins_size_overlap` is minimum insertion size overlap. Insertion size overlap is computed as:
    /// `min(len(seq1),len(seq2)) / max(len(seq1),len(seq2))`.
    ///
    /// * `min_del_size_overlap` is minimum deletion size overlap. Deletion size overlap is computed as:
    /// `min(deletion_size,deletion_size) / max(deletion_size,deletion_size)`.
    ///
    /// # Returns
    ///
    /// * Tuple of 3 VariantsList: (shared, specific to A, specific to B).
    pub fn compare(
        a: &VariantsList,
        b: &VariantsList,
        num_threads: usize,
        max_neighbor_distance: isize,
        match_all_breakpoints: bool,
        match_variant_types: bool,
        min_ins_size_overlap: f64,
        min_del_size_overlap: f64,
        variant_types_map: &HashMap<&str, String>
    ) -> (VariantsList, VariantsList, VariantsList) {
        // Step 1. Identify intersecting variants
        let vl_shared = VariantsList::intersect(
            &[a, b],
            num_threads,
            max_neighbor_distance,
            match_all_breakpoints,
            match_variant_types,
            min_ins_size_overlap,
            min_del_size_overlap,
            variant_types_map
        );

        // Step 2: Collect shared variant call IDs
        let shared_variant_call_ids: HashSet<String> = vl_shared
            .variant_call_ids()
            .into_iter()
            .collect();

        // Step 3: Helper closure to filter variants by whether they are specific to a list
        let filter_variants = |list: &VariantsList, other_set: &HashSet<String>| -> VariantsList {
            let mut filtered_list = VariantsList::new();  // Create a new empty VariantsList
            for variant in list.variants.iter() {
                if variant.variant_calls.iter().all(|vc| !other_set.contains(vc.id.as_str())) {
                    filtered_list.add_variant(variant.clone());  // Add the filtered variant to the list
                }
            }
            filtered_list
        };

        // Step 4: Identify variants specific to each list
        let vl_a_only = filter_variants(a, &shared_variant_call_ids);
        let vl_b_only = filter_variants(b, &shared_variant_call_ids);

        (vl_shared, vl_a_only, vl_b_only)
    }

    /// Filter self.variants.
    pub fn filter(&mut self, variant_filters: Vec<VariantFilter>, num_threads: usize) {
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

        // Step 3. Filter variants
        self.variants.retain(|variant| variant_ids_to_keep_set.contains(&variant.id));
    }

    /// Check if two VariantCall objects can be clustered.
    ///
    /// # Arguments
    ///
    /// * If `match_all_breakpoints==true`, both pairs of breakpoints of two variant calls
    /// must be near each other `(start1 == start2 && end1 == end2)`.
    ///
    /// * If `match_all_breakpoints==false`, only one pair of breakpoints of two variant calls
    /// must be near each other `(start1 == start2 || end1 == end2)`.
    ///
    /// * If `match_variant_types==true`, variant types (super types) must match.
    ///
    /// * `min_ins_size_overlap` is minimum insertion size overlap. Insertion size overlap is computed as:
    /// `min(len(seq1),len(seq2)) / max(len(seq1),len(seq2))`.
    ///
    /// * `min_del_size_overlap` is minimum deletion size overlap. Deletion size overlap is computed as:
    /// `min(deletion_size,deletion_size) / max(deletion_size,deletion_size)`.
    ///
    /// * `variant_types_map` is a `HashMap` where `key` is `variant type` and
    /// `value` is a `super variant type`.
    pub fn is_clusterable(
        variant_call_1: &VariantCall,
        variant_call_2: &VariantCall,
        max_neighbor_distance: isize,
        match_all_breakpoints: bool,
        match_variant_types: bool,
        min_ins_size_overlap: f64,
        min_del_size_overlap: f64,
        variant_types_map: &HashMap<&str, String>
    ) -> bool {
        let variants_list_id_1 = variant_call_1.attributes.get("temp_variants_list_id").unwrap();
        let variants_list_id_2 = variant_call_2.attributes.get("temp_variants_list_id").unwrap();
        let variant_id_1 = variant_call_1.attributes.get("temp_variant_id").unwrap();
        let variant_id_2 = variant_call_2.attributes.get("temp_variant_id").unwrap();

        // The VariantCall objects are from the same list
        // and were originally in the same Variant
        if variants_list_id_1 == variants_list_id_2 {
            if variant_id_1 == variant_id_2 && variant_call_1.id != variant_call_2.id {
                return true;
            } else {
                return false;
            }
        }

        // The VariantCall objects are from different lists
        if match_variant_types == true {
            let variant_type_1: &str = variant_types_map.get(variant_call_1.variant_type.as_str()).unwrap().as_str();
            let variant_type_2: &str = variant_types_map.get(variant_call_2.variant_type.as_str()).unwrap().as_str();

            if variant_type_1 != variant_type_2 {
                return false;
            }
            if variant_type_1 == variant_types_map[constants::DELETION].as_str() {
                let max_size: isize = max(variant_call_1.variant_size, variant_call_2.variant_size);
                let min_size: isize = min(variant_call_1.variant_size, variant_call_2.variant_size);
                let frac_size: f64 = min_size as f64 / max_size as f64;
                if frac_size < min_del_size_overlap {
                    return false;
                }
            }
            if variant_type_1 == variant_types_map[constants::INSERTION].as_str() {
                let max_size: isize = max(variant_call_1.variant_size, variant_call_2.variant_size);
                let min_size: isize = min(variant_call_1.variant_size, variant_call_2.variant_size);
                let frac_size: f64 = min_size as f64 / max_size as f64;
                if frac_size < min_ins_size_overlap {
                    return false;
                }
            }
        }

        let max_neighbor_distance_: isize;
        if variant_call_1.variant_type == constants::SINGLE_NUCLEOTIDE_VARIANT ||
            variant_call_1.variant_type == constants::MULTI_NUCLEOTIDE_VARIANT ||
            variant_call_2.variant_type == constants::SINGLE_NUCLEOTIDE_VARIANT ||
            variant_call_2.variant_type == constants::MULTI_NUCLEOTIDE_VARIANT {
            // Maximum neighbor distance must be 0 for SNVs and MNVs
            max_neighbor_distance_ = 0;
        } else {
            // Otherwise, take the input maximum neighbor distance value
            max_neighbor_distance_ = max_neighbor_distance;
        }
        let distance_11: isize = (variant_call_1.position_1 - variant_call_2.position_1).abs();
        let distance_22: isize = (variant_call_1.position_2 - variant_call_2.position_2).abs();
        let distance_12: isize = (variant_call_1.position_1 - variant_call_2.position_2).abs();
        let distance_21: isize = (variant_call_1.position_2 - variant_call_2.position_1).abs();
        if match_all_breakpoints == true {
            if ((variant_call_1.chromosome_1 == variant_call_2.chromosome_1) &&
                (variant_call_1.chromosome_2 == variant_call_2.chromosome_2) &&
                (distance_11 <= max_neighbor_distance_) &&
                (distance_22 <= max_neighbor_distance_)) ||
                ((variant_call_1.chromosome_1 == variant_call_2.chromosome_2) &&
                    (variant_call_1.chromosome_2 == variant_call_2.chromosome_1) &&
                    (distance_12 <= max_neighbor_distance_) &&
                    (distance_21 <= max_neighbor_distance_)) {
                true
            } else {
                false
            }
        } else {
            if ((variant_call_1.chromosome_1 == variant_call_2.chromosome_1) && (distance_11 <= max_neighbor_distance_)) ||
                ((variant_call_1.chromosome_2 == variant_call_2.chromosome_2) && (distance_22 <= max_neighbor_distance_)) ||
                ((variant_call_1.chromosome_1 == variant_call_2.chromosome_2) && (distance_12 <= max_neighbor_distance_)) ||
                ((variant_call_1.chromosome_2 == variant_call_2.chromosome_1) && (distance_21 <= max_neighbor_distance_)) {
                true
            } else {
                false
            }
        }
    }

    /// Identify nearby VariantCall objects.
    ///
    /// # Arguments
    ///
    /// * `variant_calls_map` is a HashMap where `key` is `
    fn identify_nearby_variant_calls(
        num_threads: usize,
        variant_calls_map: &HashMap<String, Vec<(isize, VariantCall)>>,
        max_neighbor_distance: isize,
        match_all_breakpoints: bool,
        match_variant_types: bool,
        min_ins_size_overlap: f64,
        min_del_size_overlap: f64,
        variant_types_map: &HashMap<&str, String>
    ) -> HashSet<(String,String)> {
        let thread_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();
        let results: HashSet<(String,String)> = thread_pool.install(|| {
            variant_calls_map.par_iter().flat_map(|(_chromosome, variant_call_positions)| {
                let mut pairs: HashSet<(String,String)> = HashSet::new(); // (VariantCall.id, VariantCall.id)
                for i in 0..variant_call_positions.len() {
                    for j in (i + 1)..variant_call_positions.len() {
                        let position_1: isize = variant_call_positions[i].0;
                        let position_2: isize = variant_call_positions[j].0;
                        let distance: isize = (position_1 - position_2).abs();

                        // Break the loop because subsequent VariantCall objects are outside
                        // the maximum neighbor distance
                        if distance > max_neighbor_distance {
                            break;
                        }

                        let cluster: bool = VariantsList::is_clusterable(
                            &variant_call_positions[i].1,
                            &variant_call_positions[j].1,
                            max_neighbor_distance,
                            match_all_breakpoints,
                            match_variant_types,
                            min_ins_size_overlap,
                            min_del_size_overlap,
                            variant_types_map
                        );

                        if cluster {
                            let variant_call_1_id: &str = variant_call_positions[i].1.id.as_str();
                            let variant_call_2_id: &str = variant_call_positions[j].1.id.as_str();
                            if (pairs.contains(&(variant_call_1_id.to_string(), variant_call_2_id.to_string())) == false) &&
                                (pairs.contains(&(variant_call_2_id.to_string(), variant_call_1_id.to_string())) == false) {
                                pairs.insert((variant_call_1_id.to_string(), variant_call_2_id.to_string()));
                            }
                        }
                    }
                }
                pairs
            })
            .collect()
        });
        results
    }

    /// Intersect a vector of VariantsList objects.
    ///
    /// # Arguments
    ///
    /// * If `match_all_breakpoints==true`, both pairs of breakpoints of two variant calls must be
    /// near each other `(start1 == start2 && end1 == end2)`.
    ///
    /// * If `match_all_breakpoints==false`, only one pair of breakpoints of two variant calls
    /// must be near each other `(start1 == start2 || end1 == end2)`.
    ///
    /// * If `match_variant_types==true`, variant types (super types) must match.
    ///
    /// * `min_ins_size_overlap` is minimum insertion size overlap. Insertion size overlap is computed as:
    /// `min(len(seq1),len(seq2)) / max(len(seq1),len(seq2))`.
    ///
    /// * `min_del_size_overlap` is minimum deletion size overlap. Deletion size overlap is computed as:
    /// `min(deletion_size,deletion_size) / max(deletion_size,deletion_size)`.
    ///
    /// * `variant_types_map` is `HashMap` where `key` is `variant type` and `value` is a `super variant type`.
    pub fn intersect(
        variants_lists: &[&VariantsList],
        num_threads: usize,
        max_neighbor_distance: isize,
        match_all_breakpoints: bool,
        match_variant_types: bool,
        min_ins_size_overlap: f64,
        min_del_size_overlap: f64,
        variant_types_map: &HashMap<&str, String>) -> VariantsList {
        // Step 1. Split variant calls by breakpoint chromosome
        let mut variant_calls_map: HashMap<String,Vec<(isize,VariantCall)>> = VariantsList::split_variant_calls(variants_lists);

        // Step 2. Sort variant_records_map by position
        VariantsList::sort_variant_calls(&mut variant_calls_map, num_threads);

        // Step 3. Identify nearby variant calls
        let pairs: HashSet<(String,String)> = VariantsList::identify_nearby_variant_calls(
            num_threads,
            &variant_calls_map,
            max_neighbor_distance,
            match_all_breakpoints,
            match_variant_types,
            min_ins_size_overlap,
            min_del_size_overlap,
            variant_types_map
        );

        // Step 4. By transitive property, identify clusters
        let clusters = find_clusters(pairs);

        // Step 5. Assign Variant ID to each VariantCall
        // key      =   VariantCall.id
        // value    =   Variant.id
        let mut variant_ids_map: HashMap<String, String> = HashMap::new();
        let mut new_variant_id = 1;
        for variant_call_ids in clusters.iter() {
            for variant_call_id in variant_call_ids.iter() {
                variant_ids_map.insert(variant_call_id.to_string(), new_variant_id.to_string());
            }
            new_variant_id += 1;
        }

        // Step 6. Merge VariantCall objects into common Variant objects
        let mut variants_map: HashMap<String, Variant> = HashMap::new();
        for variants_list in variants_lists.iter() {
            for variant in variants_list.variants.iter() {
                for variant_call in variant.variant_calls.iter() {
                    if let Some(variant_id) = variant_ids_map.get(&variant_call.id) {
                        // Fetch the corresponding Variant or create a new one
                        variants_map
                            .entry(variant_id.clone())
                            .or_insert_with(|| Variant::new(variant_id.clone()))
                            .add_variant_call(variant_call.clone());
                    }
                }
            }
        }

        // Step 7. Create an intersecting VariantsList object
        let mut intersecting_variants_list = VariantsList::new();
        for variant in variants_map.into_values() {
            intersecting_variants_list.add_variant(variant);
        }

        intersecting_variants_list
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
    /// * `min_ins_size_overlap`        -   Minimum insertion size overlap. Insertion size overlap is computed as:
    ///                                     min(len(seq1),len(seq2)) / max(len(seq1),len(seq2))
    /// * `variant_types_map`           -   HashMap where key is variant type and
    ///                                     value is a super set of the variant type.
    ///
    /// # Returns
    /// * `merged_variants_list`        -   Merged VariantsList object.
    pub fn merge(
        variants_lists: &[&VariantsList],
        num_threads: usize,
        max_neighbor_distance: isize,
        match_all_breakpoints: bool,
        match_variant_types: bool,
        min_ins_size_overlap: f64,
        min_del_size_overlap: f64,
        variant_types_map: &HashMap<&str, String>) -> VariantsList {
        // Step 1. Split variant calls by breakpoint chromosome
        let mut variant_calls_map: HashMap<String,Vec<(isize,VariantCall)>> = VariantsList::split_variant_calls(variants_lists);

        // Step 2. Sort variant_records_map by position
        VariantsList::sort_variant_calls(&mut variant_calls_map, num_threads);

        // Step 3. Identify nearby variant calls
        let pairs: HashSet<(String,String)> = VariantsList::identify_nearby_variant_calls(
            num_threads,
            &variant_calls_map,
            max_neighbor_distance,
            match_all_breakpoints,
            match_variant_types,
            min_ins_size_overlap,
            min_del_size_overlap,
            variant_types_map
        );

        // Step 4. By transitive property, identify clusters
        let clusters = find_clusters(pairs);

        // Step 5. Assign Variant ID to each VariantCall
        // key      =   VariantCall.id
        // value    =   Variant.id
        let mut variant_ids_map: HashMap<String, String> = HashMap::new();
        let mut new_variant_id = 1;
        for variant_call_ids in clusters.iter() {
            for variant_call_id in variant_call_ids.iter() {
                variant_ids_map.insert(variant_call_id.to_string(), new_variant_id.to_string());
            }
            new_variant_id += 1;
        }

        // Step 6. Merge VariantCall objects into common Variant objects
        let mut variants_map: HashMap<String, Variant> = HashMap::new();
        for variants_list in variants_lists.iter() {
            for variant in variants_list.variants.iter() {
                for variant_call in variant.variant_calls.iter() {
                    if let Some(variant_id) = variant_ids_map.get(&variant_call.id) {
                        // Fetch the corresponding Variant or create a new one
                        variants_map
                            .entry(variant_id.clone())
                            .or_insert_with(|| Variant::new(variant_id.clone()))
                            .add_variant_call(variant_call.clone());
                    } else {
                        // Current VariantCall does not intersect with any other VariantCall objects
                        let mut new_variant = Variant::new(new_variant_id.to_string());
                        new_variant.add_variant_call(variant_call.clone());
                        variants_map.insert(new_variant.id.clone(), new_variant);
                        new_variant_id += 1;
                    }
                }
            }
        }

        // Step 7. Create a merged VariantsList object
        let mut variants_list_merged = VariantsList::new();
        for variant in variants_map.into_values() {
            variants_list_merged.add_variant(variant);
        }

        variants_list_merged
    }

    /// Find variants overlapping query regions.
    ///
    /// # Returns
    ///
    /// * HashMap where `key` is `variant call ID` and `value` is `a vector of GenomicRange IDs`.
    pub fn overlap(
        &self,
        genomic_regions_list: GenomicRangesList,
        num_threads: usize
    ) -> HashMap<String, Vec<String>> {
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
                                let start_1 = variant_call.position_1;
                                let end_1 = variant_call.position_1;
                                let start_2 = variant_call.position_2;
                                let end_2 = variant_call.position_2;
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

        overlap_regions_map
    }

    /// Split VariantCall objects in a list of variants_list into chromsoomes.
    ///
    /// # Returns
    ///
    /// * HashMap where `key` is chromosome and `value` is `a vector of tuples (position, VariantCall)`.
    /// * Note that each VariantCall object will have a variants list ID and variant ID assigned
    /// with keys `temp_variants_list_id` and `temp_variant_id`.
    fn split_variant_calls(variants_lists: &[&VariantsList]) -> HashMap<String, Vec<(isize, VariantCall)>> {
        let mut variant_calls_map: HashMap<String, Vec<(isize, VariantCall)>> = HashMap::new();
        let mut variants_list_id: isize = 0;
        for variants_list in variants_lists.iter() {
            for variant in variants_list.variants.iter() {
                for variant_call in variant.variant_calls.iter() {
                    let mut variant_call_ = variant_call.clone();
                    variant_call_.add_attribute("temp_variants_list_id".to_string(), variants_list_id.to_string());
                    variant_call_.add_attribute("temp_variant_id".to_string(), variant.id.to_string());

                    // Insert position 1
                    variant_calls_map
                        .entry(variant_call_.chromosome_1.clone())
                        .or_insert_with(Vec::new)
                        .push((variant_call_.position_1, variant_call_.clone()));

                    // Insert position 2 (if on different chromosomes)
                    variant_calls_map
                        .entry(variant_call_.chromosome_2.clone())
                        .or_insert_with(Vec::new)
                        .push((variant_call_.position_2, variant_call_.clone()));
                }
            }
            variants_list_id += 1;
        }
        variant_calls_map
    }

    /// Sorts `self.variants` by the first `VariantCall` object's `chromosome_1` and `position_1`.
    pub fn sort(&mut self) {
        // Sorting by the first VariantCall's chromosome_1 and position_1
        self.variants.sort_by(|a, b| {
            let a_first_call = &a.variant_calls[0];
            let b_first_call = &b.variant_calls[0];

            // First, compare by chromosome_1
            let cmp = a_first_call.chromosome_1.cmp(&b_first_call.chromosome_1);

            // If chromosomes are the same, compare by position_1
            if cmp == std::cmp::Ordering::Equal {
                a_first_call.position_1.cmp(&b_first_call.position_1)
            } else {
                cmp
            }
        });
    }

    /// Sorts variant calls by position.
    ///
    /// # Arguments
    ///
    /// * `variant_calls_map` is a HashMap where the `key` is `chromosome` and the `value` is
    /// `a vector of tuples (position, VariantCall)`.
    fn sort_variant_calls(
        variant_calls_map: &mut HashMap<String, Vec<(isize, VariantCall)>>,
        num_threads: usize
    ) {
        let thread_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();
        thread_pool.install(|| {
            variant_calls_map.par_iter_mut().for_each(|(_key, variant_calls)| {
                variant_calls.sort_by_key(|&(position, _)| position);
            });
        });
    }

    /// Subtract another VariantsList object.
    ///
    /// # Arguments
    ///
    pub fn subtract(
        &self,
        variants_list: &VariantsList,
        num_threads: usize,
        max_neighbor_distance: isize,
        match_all_breakpoints: bool,
        match_variant_types: bool,
        min_ins_size_overlap: f64,
        min_del_size_overlap: f64,
        variant_types_map: &HashMap<&str, String>
    ) -> VariantsList {
        let (vl_shared, vl_a_only, vl_b_only) = VariantsList::compare(
            &self,
            variants_list,
            num_threads,
            max_neighbor_distance,
            match_all_breakpoints,
            match_variant_types,
            min_ins_size_overlap,
            min_del_size_overlap,
            variant_types_map
        );
        vl_a_only
    }

    pub fn variant_call_ids(&self) -> Vec<String> {
        let mut variant_call_ids: Vec<String> = Vec::new();
        for variant in self.variants.iter() {
            for variant_call in variant.variant_calls.iter() {
                variant_call_ids.push(variant_call.id.to_string());
            }
        }
        variant_call_ids
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
