# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


"""
The purpose of this python3 script is to implement the PositionAnnotation dataclass.
"""


import pandas as pd
from dataclasses import dataclass, field


@dataclass
class PositionAnnotation:
    annotator: str
    region: str
    species: str
    annotator_version: str = field(default='')
    gene_id: str = field(default='')
    gene_id_stable: str = field(default='')
    gene_name: str = field(default='')
    gene_strand: str = field(default='')
    gene_type: str = field(default='')
    gene_version: str = field(default='')
    transcript_id: str = field(default='')
    transcript_id_stable: str = field(default='')
    transcript_name: str = field(default='')
    transcript_strand: str = field(default='')
    transcript_type: str = field(default='')
    transcript_version: str = field(default='')
    exon_id: str = field(default='')
    exon_id_stable: str = field(default='')

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self.to_dataframe_row())

    def to_dataframe_row(self):
        return {
            'annotator': [self.annotator],
            'annotator_version': [self.annotator_version],
            'gene_id': [self.gene_id],
            'gene_id_stable': [self.gene_id_stable],
            'gene_name': [self.gene_name],
            'gene_strand': [self.gene_strand],
            'gene_type': [self.gene_type],
            'gene_version': [self.gene_version],
            'transcript_id': [self.transcript_id],
            'transcript_id_stable': [self.transcript_id_stable],
            'transcript_name': [self.transcript_name],
            'transcript_strand': [self.transcript_strand],
            'transcript_type': [self.transcript_type],
            'transcript_version': [self.transcript_version],
            'exon_id': [self.exon_id],
            'exon_id_stable': [self.exon_id_stable],
            'region': [self.region],
            'species': [self.species]
        }

    def to_dict(self):
        return {
            'annotator': self.annotator,
            'annotator_version': self.annotator_version,
            'gene_id': self.gene_id,
            'gene_id_stable': self.gene_id_stable,
            'gene_name': self.gene_name,
            'gene_strand': self.gene_strand,
            'gene_type': self.gene_type,
            'gene_version': self.gene_version,
            'transcript_id': self.transcript_id,
            'transcript_id_stable': self.transcript_id_stable,
            'transcript_name': self.transcript_name,
            'transcript_strand': self.transcript_strand,
            'transcript_type': self.transcript_type,
            'transcript_version': self.transcript_version,
            'exon_id': self.exon_id,
            'exon_id_stable': self.exon_id_stable,
            'region': self.region,
            'species': self.species
        }
