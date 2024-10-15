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


extern crate chrono;
extern crate env_logger;
extern crate pyo3;
extern crate serde_json;
extern crate vstol;
use chrono::Local;
use env_logger::{Builder, Env};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use std::collections::HashMap;
use std::io::Write;
use vstol::constants;
use vstol::metrics::calculate_average_alignment_scores as calculate_average_alignment_scores_;
use vstol::structs::genomic_range::GenomicRange;
use vstol::structs::genomic_ranges_list::GenomicRangesList;
use vstol::structs::variant::Variant;
use vstol::structs::variant_call::VariantCall;
use vstol::structs::variant_call_annotation::VariantCallAnnotation;
use vstol::structs::variant_filter::VariantFilter;
use vstol::structs::variants_list::VariantsList;

mod interface;
use crate::interface::*;


/// Calculates average alignment scores for a given list of regions.
///
/// # Arguments
/// * `bam_file`        -   BAM file.
/// * `regions`         -   vector of tuples where each tuple is (chromosome,start,end).
/// * `num_threads`     -   number of threads.
///
/// # Returns
/// * HashMap where key = (chromosome,start,end) and value = average alignment score.
#[pyfunction]
fn calculate_average_alignment_scores(
    py: Python,
    bam_file: String,
    regions: Vec<(String,u32,u32)>,
    num_threads: usize) -> PyResult<HashMap<(String, u32, u32), f64>> {
    let scores: HashMap<(String,u32,u32),f64> = calculate_average_alignment_scores_(
        bam_file.as_str(),
        &regions,
        num_threads
    );
    Ok(scores)
}

/// This function compares two serialized VariantsList objects and returns three VariantsList objects.
///
/// # Arguments
/// * `py_list`                 -   list of 2 serialized VariantsList objects.
/// * `num_threads`             -   number of threads.
///
/// # Returns
/// * A serialized list of 3 VariantsList strings.
#[pyfunction]
fn compare_variants_lists(
    py: Python,
    vl_list: &PyList,
    num_threads: usize,
    max_neighbor_distance: isize,
    match_all_breakpoints: bool,
    match_variant_types: bool,
    min_ins_size_overlap: f64,
    min_del_size_overlap: f64
) -> PyResult<String> {
    // Step 1. Deserialize VariantsList objects
    let variants_lists: Vec<VariantsList> = deserialize_variants_lists(&vl_list);

    assert_eq!(variants_lists.len(), 2);

    // Step 2. Compare VariantsList objects
    let (vl_shared, vl_a_only, vl_b_only) = VariantsList::compare(
        &variants_lists[0],
        &variants_lists[1],
        num_threads,
        max_neighbor_distance,
        match_all_breakpoints,
        match_variant_types,
        min_ins_size_overlap,
        min_del_size_overlap,
        &constants::VARIANT_TYPES_MAP
    );

    // Step 3. Serialize VariantsList objects
    let serialized = serde_json::to_string(&[vl_shared, vl_a_only, vl_b_only]).expect("Serialization of the list of VariantsList objects failed");

    Ok(serialized)
}

/// This function filters a serialized VariantsList object and returns a filtered VariantsList.
///
/// # Arguments
/// * `vl_target`               -   serialized VariantsList object.
/// * `vl_target`               -   list of serialized VariantFilter objects.
/// * `num_threads`             -   number of threads.
///
/// # Returns
/// * A serialized VariantsList string.
#[pyfunction]
fn filter_variants_list(
    py: Python,
    vl_target: String,
    filter_list: &PyList,
    num_threads: usize
) -> PyResult<String> {
    // Step 1. Deserialize VariantsList object
    let mut variants_list: VariantsList = deserialize_variants_list(&vl_target);

    // Step 2. Deserialize VariantFilter objects
    let variant_filters: Vec<VariantFilter> = deserialize_variant_filters(filter_list);

    // Step 3. Filter VariantsList object
    variants_list.filter(variant_filters, num_threads);

    // Step 4. Serialize filtered VariantsList object
    let serialized = serde_json::to_string(&variants_list).expect("Serialization of the filtered VariantsList object failed");

    Ok(serialized)
}

/// This function identifies intersecting (or nearby) variant calls given
/// a vector of serialized VariantsList objects.
///
/// # Arguments
/// * `vl_list`                 -   list of serialized VariantsList objects.
/// * `num_threads`             -   number of threads.
/// * `max_neighbor_distance`   -   maximum neighbor distance.
/// * `match_all_breakpoints`   -   If true, both pairs of breakpoints of two variant calls
///                                 must be near each other (start1 == start2 && end1 == end2).
///                                 If false, only one pair of breakpoints of two variant calls
///                                 must be near each other (start1 == start2 || end1 == end2).
/// * `match_variant_types`     -   If true, variant types (super types) must match.
///
/// # Returns
/// * A serialized VariantsList object.
#[pyfunction]
fn intersect_variants_lists(
    py: Python,
    vl_list: &PyList,
    num_threads: usize,
    max_neighbor_distance: isize,
    match_all_breakpoints: bool,
    match_variant_types: bool,
    min_ins_size_overlap: f64,
    min_del_size_overlap: f64
) -> PyResult<String> {
    // Step 1. Deserialize VariantsList objects
    let variants_lists: Vec<VariantsList> = deserialize_variants_lists(vl_list);
    let variants_refs: Vec<&VariantsList> = variants_lists.iter().collect();

    // Step 2. Identify intersecting (or nearby) variant calls in VariantsList objects
    let intersecting_variants_list: VariantsList = VariantsList::intersect(
        &variants_refs,
        num_threads,
        max_neighbor_distance,
        match_all_breakpoints,
        match_variant_types,
        min_ins_size_overlap,
        min_del_size_overlap,
        &constants::VARIANT_TYPES_MAP
    );

    // Step 3. Serialize merged VariantsList object
    let serialized = serde_json::to_string(&intersecting_variants_list).expect("Serialization of the merged VariantsList object failed");

    Ok(serialized)
}

/// This function merges a vector of serialized VariantsList objects into one.
///
/// # Arguments
/// * `vl_list`                 -   list of serialized VariantsList objects.
/// * `num_threads`             -   number of threads.
/// * `max_neighbor_distance`   -   maximum neighbor distance.
///
/// # Returns
/// * A serialized VariantsList object.
#[pyfunction]
fn merge_variants_lists(
    py: Python,
    vl_list: &PyList,
    num_threads: usize,
    max_neighbor_distance: isize,
    match_all_breakpoints: bool,
    match_variant_types: bool,
    min_ins_size_overlap: f64,
    min_del_size_overlap: f64) -> PyResult<String> {
    // Step 1. Deserialize VariantsList objects
    let variants_lists: Vec<VariantsList> = deserialize_variants_lists(vl_list);
    let variants_refs: Vec<&VariantsList> = variants_lists.iter().collect();

    // Step 2. Merge VariantsList objects
    let merged_variants_list: VariantsList = VariantsList::merge(
        &variants_refs,
        num_threads,
        max_neighbor_distance,
        match_all_breakpoints,
        match_variant_types,
        min_ins_size_overlap,
        min_del_size_overlap,
        &constants::VARIANT_TYPES_MAP
    );

    // Step 3. Serialize merged VariantsList object
    let serialized = serde_json::to_string(&merged_variants_list).expect("Serialization of the merged VariantsList object failed");

    Ok(serialized)
}

/// This function identifies overlapping VariantCall objects.
///
/// # Arguments
/// * `vl_target`                       -   serialized VariantsList object.
/// * `granges_list`                    -   serialized GenomicRangesList object.
/// * `num_threads`                     -   number of threads.
///
/// # Returns
/// * `overlapping_variant_call_ids`    -   HashMap where key is VariantCall.id and
///                                         value is a vector of GenomicRange.id
#[pyfunction]
fn overlap_variant_calls(
    py: Python,
    vl_target: String,
    granges_list: String,
    num_threads: usize
) -> Py<PyAny> {
    // Step 1. Deserialize VariantsList object
    let mut variants_list: VariantsList = deserialize_variants_list(&vl_target);

    // Step 2. Deserialize GenomicRangesList objects
    let genomic_ranges_list: GenomicRangesList = deserialize_genomic_ranges_list(&granges_list);

    // Step 3. Find overlapping VariantCall IDs
    let overlapping_variant_call_ids: HashMap<String, Vec<String>> = variants_list.overlap(
        genomic_ranges_list,
        num_threads
    );

    return Python::with_gil(|py| {
        overlapping_variant_call_ids.to_object(py)
    });
}

/// This function merges a vector of serialized VariantsList objects into one.
///
/// # Arguments
///
/// * `a`                       -   serialized VariantsList object.
/// * `b`                       -   serialized VariantsList object.
/// * `num_threads`             -   number of threads.
/// * `max_neighbor_distance`   -   maximum neighbor distance.
/// * `match_variant_types`     -   If true, match variant types.
/// * `match_all_breakpoints`   -   If true, match all breakpoints.
/// * `min_ins_size_overlap`    -   Minimum insertion size overlap.
/// * `min_del_size_overlap`    -   Minimum deletion size overlap.
///
/// # Returns
///
/// * A serialized VariantsList object.
#[pyfunction]
fn subtract_variants_list(
    py: Python,
    a: String,
    b: String,
    num_threads: usize,
    max_neighbor_distance: isize,
    match_all_breakpoints: bool,
    match_variant_types: bool,
    min_ins_size_overlap: f64,
    min_del_size_overlap: f64) -> PyResult<String> {
    // Step 1. Deserialize VariantsList objects
    let vl_a: VariantsList = deserialize_variants_list(&a);
    let vl_b: VariantsList = deserialize_variants_list(&b);

    // Step 2. Subtract VariantsList object B from A
    let vl_subtracted: VariantsList = vl_a.subtract(
        &vl_b,
        num_threads,
        max_neighbor_distance,
        match_all_breakpoints,
        match_variant_types,
        min_ins_size_overlap,
        min_del_size_overlap,
        &constants::VARIANT_TYPES_MAP
    );

    // Step 3. Serialize subtracted VariantsList object
    let serialized = serde_json::to_string(&vl_subtracted).expect("Serialization of the subtracted VariantsList object failed");

    Ok(serialized)
}

#[pymodule]
fn vstolibrs(_py: Python, m: &PyModule) -> PyResult<()> {
    // Initialize the logger
    Builder::from_env(Env::default().default_filter_or("info")).format(|buf, record| {
        let now = Local::now();
        writeln!(
            buf,
            "{} {} [{:>50}] {}",
            now.format("%Y-%m-%d %H:%M:%S"),
            record.level(),
            record.target(),
            record.args()
        )
    }).init();

    m.add_function(wrap_pyfunction!(calculate_average_alignment_scores, m)?);
    m.add_function(wrap_pyfunction!(compare_variants_lists, m)?);
    m.add_function(wrap_pyfunction!(filter_variants_list, m)?);
    m.add_function(wrap_pyfunction!(intersect_variants_lists, m)?);
    m.add_function(wrap_pyfunction!(merge_variants_lists, m)?);
    m.add_function(wrap_pyfunction!(overlap_variant_calls, m)?);
    m.add_function(wrap_pyfunction!(subtract_variants_list, m)?);
    Ok(())
}
