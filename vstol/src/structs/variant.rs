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


extern crate exitcode;
extern crate serde;
use serde::{Serialize, Deserialize};
use crate::common::{calculate_max,calculate_mean,calculate_median,calculate_min};
use crate::structs::variant_call::VariantCall;


#[derive(Debug, Serialize, Deserialize)]
pub struct Variant {
    pub id: String,
    pub variant_calls: Vec<VariantCall>
}

impl Variant {
    pub fn new(id: String) -> Self {
        Self {
            id: id,
            variant_calls: Vec::new()
        }
    }

    pub fn add_variant_call(&mut self, variant_call: VariantCall) {
        // Find an index position to insert the new VariantCall object
        let insert_idx = self.variant_calls.binary_search_by(|item| item.position_1.cmp(&variant_call.position_1));
        match insert_idx {
            Ok(idx) => {
                self.variant_calls.insert(idx, variant_call);
            }
            Err(idx) => {
                self.variant_calls.insert(idx, variant_call);
            }
        }
    }

    pub fn get_mean(&self, attribute: String) -> f64 {
        let values: Vec<f64> = self.variant_calls.iter().map(|variant_call| {
            let value: f64 = match attribute.as_str() {
                "position_1" => variant_call.position_1 as f64,
                "position_2" => variant_call.position_2 as f64,
                "quality_score" => variant_call.quality_score as f64,
                "variant_size" => variant_call.variant_size as f64,
                "total_read_count" => variant_call.total_read_count as f64,
                "reference_allele_read_count" => variant_call.reference_allele_read_count as f64,
                "alternate_allele_read_count" => variant_call.alternate_allele_read_count as f64,
                "alternate_allele_fraction" => variant_call.alternate_allele_fraction as f64,
                _ => -1.0,
            };
            return value;
        }).collect();
        calculate_mean(values)
    }

    pub fn get_median(&self, attribute: String) -> f64 {
        let values: Vec<f64> = self.variant_calls.iter().map(|variant_call| {
            let value: f64 = match attribute.as_str() {
                "position_1" => variant_call.position_1 as f64,
                "position_2" => variant_call.position_2 as f64,
                "quality_score" => variant_call.quality_score as f64,
                "variant_size" => variant_call.variant_size as f64,
                "total_read_count" => variant_call.total_read_count as f64,
                "reference_allele_read_count" => variant_call.reference_allele_read_count as f64,
                "alternate_allele_read_count" => variant_call.alternate_allele_read_count as f64,
                "alternate_allele_fraction" => variant_call.alternate_allele_fraction as f64,
                _ => -1.0,
            };
            return value;
        }).collect();
        calculate_median(values)
    }

    pub fn get_max(&self, attribute: String) -> f64 {
        let values: Vec<f64> = self.variant_calls.iter().map(|variant_call| {
            let value: f64 = match attribute.as_str() {
                "position_1" => variant_call.position_1 as f64,
                "position_2" => variant_call.position_2 as f64,
                "quality_score" => variant_call.quality_score as f64,
                "variant_size" => variant_call.variant_size as f64,
                "total_read_count" => variant_call.total_read_count as f64,
                "reference_allele_read_count" => variant_call.reference_allele_read_count as f64,
                "alternate_allele_read_count" => variant_call.alternate_allele_read_count as f64,
                "alternate_allele_fraction" => variant_call.alternate_allele_fraction as f64,
                _ => -1.0,
            };
            return value;
        }).collect();
        calculate_max(values)
    }

    pub fn get_min(&self, attribute: String) -> f64 {
        let values: Vec<f64> = self.variant_calls.iter().map(|variant_call| {
            let value: f64 = match attribute.as_str() {
                "position_1" => variant_call.position_1 as f64,
                "position_2" => variant_call.position_2 as f64,
                "quality_score" => variant_call.quality_score as f64,
                "variant_size" => variant_call.variant_size as f64,
                "total_read_count" => variant_call.total_read_count as f64,
                "reference_allele_read_count" => variant_call.reference_allele_read_count as f64,
                "alternate_allele_read_count" => variant_call.alternate_allele_read_count as f64,
                "alternate_allele_fraction" => variant_call.alternate_allele_fraction as f64,
                _ => -1.0,
            };
            return value;
        }).collect();
        calculate_min(values)
    }
}

impl Clone for Variant {
    fn clone(&self) -> Self {
        let mut variant_calls: Vec<VariantCall> = Vec::new();
        for variant_call in &self.variant_calls {
            variant_calls.push(variant_call.clone());
        }
        Variant {
            id: self.id.clone(),
            variant_calls: variant_calls
        }
    }
}