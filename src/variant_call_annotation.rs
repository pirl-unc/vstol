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


extern crate serde;
use serde::{Serialize, Deserialize};


#[derive(Debug, Serialize, Deserialize)]
pub struct VariantCallAnnotation {
    annotator: String,
    annotator_version: String,
    gene_id: String,
    gene_id_stable: String,
    gene_name: String,
    gene_strand: String,
    gene_type: String,
    gene_version: String,
    region: String,
    species: String
}

impl VariantCallAnnotation {
    fn new(
        annotator: String,
        annotator_version: String,
        gene_id: String,
        gene_id_stable: String,
        gene_name: String,
        gene_strand: String,
        gene_type: String,
        gene_version: String,
        region: String,
        species: String) -> Self {
        VariantCallAnnotation {
            annotator: annotator,
            annotator_version: annotator_version,
            gene_id: gene_id,
            gene_id_stable: gene_id_stable,
            gene_name: gene_name,
            gene_strand: gene_strand,
            gene_type: gene_type,
            gene_version: gene_version,
            region: region,
            species: species
        }
    }
}

impl Clone for VariantCallAnnotation {
    fn clone(&self) -> Self {
        VariantCallAnnotation {
            annotator: self.annotator.clone(),
            annotator_version: self.annotator_version.clone(),
            gene_id: self.gene_id.clone(),
            gene_id_stable: self.gene_id_stable.clone(),
            gene_name: self.gene_name.clone(),
            gene_strand: self.gene_strand.clone(),
            gene_type: self.gene_type.clone(),
            gene_version: self.gene_version.clone(),
            region: self.region.clone(),
            species: self.species.clone()
        }
    }
}
