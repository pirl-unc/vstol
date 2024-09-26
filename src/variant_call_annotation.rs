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
    transcript_id: String,
    transcript_id_stable: String,
    transcript_name: String,
    transcript_strand: String,
    transcript_type: String,
    transcript_version: String,
    exon_id: String,
    exon_id_stable: String,
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
        transcript_id: String,
        transcript_id_stable: String,
        transcript_name: String,
        transcript_strand: String,
        transcript_type: String,
        transcript_version: String,
        exon_id: String,
        exon_id_stable: String,
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
            transcript_id: transcript_id,
            transcript_id_stable: transcript_id_stable,
            transcript_name: transcript_name,
            transcript_strand: transcript_strand,
            transcript_type: transcript_type,
            transcript_version: transcript_version,
            exon_id: exon_id,
            exon_id_stable: exon_id_stable,
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
            transcript_id: self.transcript_id.clone(),
            transcript_id_stable: self.transcript_id_stable.clone(),
            transcript_name: self.transcript_name.clone(),
            transcript_strand: self.transcript_strand.clone(),
            transcript_type: self.transcript_type.clone(),
            transcript_version: self.transcript_version.clone(),
            exon_id: self.exon_id.clone(),
            exon_id_stable: self.exon_id_stable.clone(),
            region: self.region.clone(),
            species: self.species.clone()
        }
    }
}
