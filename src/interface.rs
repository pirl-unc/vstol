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


extern crate pyo3;
extern crate serde_json;
extern crate vstol;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use std::collections::HashMap;
use vstol::structs::genomic_range::GenomicRange;
use vstol::structs::genomic_ranges_list::GenomicRangesList;
use vstol::structs::variant::Variant;
use vstol::structs::variant_call::VariantCall;
use vstol::structs::variant_call_annotation::VariantCallAnnotation;
use vstol::structs::variant_filter::VariantFilter;
use vstol::structs::variants_list::VariantsList;


/// This function deserializes a serialized GenomicRangesList object.
///
/// # Arguments
/// * `json_str`                        -   serialized GenomicRangesList object.
///
/// # Returns
/// * `deserialize_genomic_ranges_list` -   GenomicRangesList object.
pub fn deserialize_genomic_ranges_list(json_str: &str) -> GenomicRangesList {
    let genomic_ranges_list: Result<GenomicRangesList, serde_json::Error> = serde_json::from_str(json_str);
    match genomic_ranges_list {
        Ok(result) => {
            return result;
        }
        Err(e) => {
            panic!("Error deserializing JSON: {}", e);
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
pub fn deserialize_variants_lists(py_list: &PyList) -> Vec<VariantsList> {
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
pub fn deserialize_variants_list(json_str: &str) -> VariantsList {
    let variants_list: Result<VariantsList, serde_json::Error> = serde_json::from_str(json_str);
    match variants_list {
        Ok(result) => {
            return result;
        }
        Err(e) => {
            panic!("Error deserializing JSON: {}", e);
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
pub fn deserialize_variant_filters(py_list: &PyList) -> Vec<VariantFilter> {
    let mut variant_filters: Vec<VariantFilter> = Vec::new();
    for py_str in py_list.iter() {
        variant_filters.push(deserialize_variant_filter(&py_str.to_string()));
    }
    return variant_filters;
}

pub fn deserialize_variant_filter(json_str: &str) -> VariantFilter {
    let variant_filter: Result<VariantFilter, serde_json::Error> = serde_json::from_str(json_str);
    match variant_filter {
        Ok(result) => {
            return result;
        }
        Err(e) => {
            panic!("Error deserializing JSON: {}", e);
        }
    }
}