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
extern crate exitcode;
extern crate log;
extern crate pyo3;
extern crate serde_json;
use chrono::Local;
use env_logger::{Builder, Env};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use std::collections::HashMap;
use std::io::Write;
mod constants;
mod genomic_range;
mod genomic_ranges_list;
mod metrics;
mod utilities;
mod variant;
mod variant_call;
mod variant_call_annotation;
mod variant_filter;
mod variants_list;
use genomic_range::GenomicRange;
use genomic_ranges_list::GenomicRangesList;
use metrics::calculate_average_alignment_scores as calculate_average_alignment_scores_;
use variant::Variant;
use variant_call::VariantCall;
use variant_call_annotation::VariantCallAnnotation;
use variant_filter::VariantFilter;
use variants_list::VariantsList;


/// This function deserializes a serialized GenomicRangesList object.
///
/// # Arguments
/// * `json_str`                        -   serialized GenomicRangesList object.
///
/// # Returns
/// * `deserialize_genomic_ranges_list` -   GenomicRangesList object.
fn deserialize_genomic_ranges_list(json_str: &str) -> GenomicRangesList {
    let genomic_ranges_list: Result<GenomicRangesList, serde_json::Error> = serde_json::from_str(json_str);
    match genomic_ranges_list {
        Ok(result) => {
            return result;
        }
        Err(e) => {
            eprintln!("Error deserializing JSON: {}", e);
            std::process::exit(exitcode::DATAERR);
        }
    }
}

/// This function deserializes a list of serialized VariantsList objects.
///
/// # Arguments
/// * `py_list`         -   list of serialized VariantsList objects.
///
/// # Returns
/// * `variants_list`   -   vector of VariantsList objects.
fn deserialize_variants_lists(py_list: &PyList) -> Vec<VariantsList> {
    let mut variants_lists: Vec<VariantsList> = Vec::new();
    for py_str in py_list.iter() {
        variants_lists.push(deserialize_variants_list(&py_str.to_string()));
    }
    return variants_lists;
}

/// This function deserializes a VariantsList object.
///
/// # Arguments
/// * `json_str`        -   serialized VariantsList object.
///
/// # Returns
/// * `variants_list`   -   VariantsList object.
fn deserialize_variants_list(json_str: &str) -> VariantsList {
    let variants_list: Result<VariantsList, serde_json::Error> = serde_json::from_str(json_str);
    match variants_list {
        Ok(result) => {
            return result;
        }
        Err(e) => {
            eprintln!("Error deserializing JSON: {}", e);
            std::process::exit(exitcode::DATAERR);
        }
    }
}

/// This function deserializes a list of VariantFilter objects.
///
/// # Arguments
/// * `py_list`             -   a list of serialized VariantsFilter objects.
///
/// # Returns
/// * `variant_filters`     -   a vector of VariantFilter objects.
fn deserialize_variant_filters(py_list: &PyList) -> Vec<VariantFilter> {
    let mut variant_filters: Vec<VariantFilter> = Vec::new();
    for py_str in py_list.iter() {
        variant_filters.push(deserialize_variant_filter(&py_str.to_string()));
    }
    return variant_filters;
}

fn deserialize_variant_filter(json_str: &str) -> VariantFilter {
    let variant_filter: Result<VariantFilter, serde_json::Error> = serde_json::from_str(json_str);
    match variant_filter {
        Ok(result) => {
            return result;
        }
        Err(e) => {
            eprintln!("Error deserializing JSON: {}", e);
            std::process::exit(exitcode::DATAERR);
        }
    }
}

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

/// This function filters a serialized VariantsList object and returns a filtered VariantsList.
///
/// # Arguments
/// * `py_str`                  -   serialized VariantsList object.
/// * `py_list`                 -   list of serialized VariantFilter objects.
/// * `num_threads`             -   number of threads.
///
/// # Returns
/// * A serialized VariantsList string.
#[pyfunction]
fn filter_variants_list(
    py: Python,
    py_str: String,
    py_list: &PyList,
    num_threads: usize) -> PyResult<String> {
    // Step 1. Deserialize VariantsList object
    let mut variants_list: VariantsList = deserialize_variants_list(&py_str);

    // Step 2. Deserialize VariantFilter objects
    let variant_filters: Vec<VariantFilter> = deserialize_variant_filters(py_list);

    // Step 3. Filter VariantsList object
    let filtered_variants_list: VariantsList = variants_list.filter(variant_filters, num_threads);

    // Step 4. Serialize filtered VariantsList object
    let serialized = serde_json::to_string(&filtered_variants_list).expect("Serialization of the filtered VariantsList object failed");

    Ok(serialized)
}

/// This function identifies overlapping VariantCall objects.
///
/// # Arguments
/// * `variants_list_str`               -   serialized VariantsList object.
/// * `genomic_ranges_list_str`         -   serialized GenomicRangesList object.
/// * `num_threads`                     -   number of threads.
/// * `padding`                         -   padding to apply to GenomicRange start and end.
///
/// # Returns
/// * `overlapping_variant_call_ids`    -   HashMap where key is VariantCall.id and
///                                         value is a vector of GenomicRange.id
#[pyfunction]
fn find_overlapping_variant_calls(
    py: Python,
    variants_list_str: String,
    genomic_ranges_list_str: String,
    num_threads: usize,
    padding: isize) -> Py<PyAny> {
    // Step 1. Deserialize VariantsList object
    let mut variants_list: VariantsList = deserialize_variants_list(&variants_list_str);

    // Step 2. Deserialize GenomicRangesList objects
    let genomic_ranges_list: GenomicRangesList = deserialize_genomic_ranges_list(&genomic_ranges_list_str);

    // Step 3. Find overlapping VariantCall IDs
    let overlapping_variant_call_ids: HashMap<String, Vec<String>> = variants_list.overlap(
        genomic_ranges_list,
        num_threads,
        padding
    );

    return Python::with_gil(|py| {
        overlapping_variant_call_ids.to_object(py)
    });
}

/// This function identifies intersecting (or nearby) variant calls given
/// a vector of serialized VariantsList objects.
///
/// # Arguments
/// * `py_list`                 -   list of serialized VariantsList objects.
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
    py: Python, py_list: &PyList,
    num_threads: usize,
    max_neighbor_distance: isize,
    match_all_breakpoints: bool,
    match_variant_types: bool) -> PyResult<String> {
    // Step 1. Deserialize VariantsList objects
    let mut variants_lists: Vec<VariantsList> = deserialize_variants_lists(py_list);

    // Step 2. Identify intersecting (or nearby) variant calls in VariantsList objects
    let intersecting_variants_list: VariantsList = VariantsList::merge(
        variants_lists,
        num_threads,
        max_neighbor_distance,
        match_all_breakpoints,
        match_variant_types,
        true,
        &constants::VARIANT_TYPES_MAP
    );

    // Step 3. Serialize merged VariantsList object
    let serialized = serde_json::to_string(&intersecting_variants_list).expect("Serialization of the merged VariantsList object failed");

    Ok(serialized)
}

/// This function merges a vector of serialized VariantsList objects into one.
///
/// # Arguments
/// * `py_list`                 -   list of serialized VariantsList objects.
/// * `num_threads`             -   number of threads.
/// * `max_neighbor_distance`   -   maximum neighbor distance.
///
/// # Returns
/// * A serialized VariantsList object.
#[pyfunction]
fn merge_variants_lists(
    py: Python, py_list: &PyList,
    num_threads: usize,
    max_neighbor_distance: isize,
    match_all_breakpoints: bool,
    match_variant_types: bool) -> PyResult<String> {
    // Step 1. Deserialize VariantsList objects
    let mut variants_lists: Vec<VariantsList> = deserialize_variants_lists(py_list);

    // Step 2. Merge VariantsList objects
    let merged_variants_list: VariantsList = VariantsList::merge(
        variants_lists,
        num_threads,
        max_neighbor_distance,
        match_all_breakpoints,
        match_variant_types,
        false,
        &constants::VARIANT_TYPES_MAP
    );

    // Step 3. Serialize merged VariantsList object
    let serialized = serde_json::to_string(&merged_variants_list).expect("Serialization of the merged VariantsList object failed");

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
    m.add_function(wrap_pyfunction!(filter_variants_list, m)?);
    m.add_function(wrap_pyfunction!(find_overlapping_variant_calls, m)?);
    m.add_function(wrap_pyfunction!(intersect_variants_lists, m)?);
    m.add_function(wrap_pyfunction!(merge_variants_lists, m)?);
    Ok(())
}
