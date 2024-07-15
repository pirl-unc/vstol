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
The purpose of this python3 script is to implement the VariantCall dataclass.
"""


import re
import pandas as pd
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Dict, List, Set, Tuple
from functools import total_ordering
from .constants import TranslocationOrientations, VariantTypes
from .variant_call_annotation import VariantCallAnnotation


@total_ordering
@dataclass
class VariantCall:
    # Mandatory fields
    id: str
    sample_id: str
    chromosome_1: str
    position_1: int
    chromosome_2: str
    position_2: int
    variant_type: str
    reference_allele: str
    alternate_allele: str

    # Optional fields
    source_id: str = field(default='')
    phase_block_id: str = field(default='')
    clone_id: str = field(default='')
    nucleic_acid: str = field(default='')
    variant_calling_method: str = field(default='')
    sequencing_platform: str = field(default='')
    filter: str = field(default='')
    quality_score: float = field(default=-1.0)
    precise: bool = field(default=False)
    variant_subtype: str = field(default='')
    variant_size: int = field(default=-1)
    reference_allele_read_count: int = field(default=-1)
    alternate_allele_read_count: int = field(default=-1)
    total_read_count: int = field(default=-1)
    alternate_allele_fraction: float = field(default=-1.0)
    alternate_allele_read_ids: Set[str] = field(default_factory=set)
    variant_sequences: Set[str] = field(default_factory=set)
    attributes: OrderedDict = field(default_factory=dict)
    tags: Set[str] = field(default_factory=set)
    average_alignment_score_window: int = field(default=-1)
    position_1_average_alignment_score: float = field(default=-1.0)
    position_2_average_alignment_score: float = field(default=-1.0)
    position_1_annotations: List[VariantCallAnnotation] = field(default_factory=list)
    position_2_annotations: List[VariantCallAnnotation] = field(default_factory=list)

    def __lt__(self, other):
        if isinstance(other, VariantCall):
            return (self.chromosome_1,
                    self.position_1,
                    self.chromosome_2,
                    self.position_2) < \
                   (other.chromosome_1,
                    other.position_1,
                    other.chromosome_2,
                    other.position_2)
        return NotImplemented

    def __eq__(self, other):
        if isinstance(other, VariantCall):
            return (self.chromosome_1,
                    self.position_1,
                    self.chromosome_2,
                    self.position_2) == \
                   (other.chromosome_1,
                    other.position_1,
                    other.chromosome_2,
                    other.position_2)
        return NotImplemented

    def add_position_1_annotation(self, variant_call_annotation: VariantCallAnnotation):
        self.position_1_annotations.append(variant_call_annotation)

    def add_position_2_annotation(self, variant_call_annotation: VariantCallAnnotation):
        self.position_2_annotations.append(variant_call_annotation)

    def get_translocation_orientation(self) -> Tuple[str,str,int,str,int]:
        """
        Get translocation orientation.

        Returns:
            Tuple[orientation,t_chromosome,t_position,p_chromosome,p_position]
        """
        if self.variant_type == VariantTypes.TRANSLOCATION:
            if re.search("^.*\[.*\[$", self.alternate_allele):                  # t[p[ piece extending to the right of p is joined after t
                orientation = TranslocationOrientations.ORIENTATION_1
                alternate_allele_elements = self.alternate_allele.split('[')
                t = alternate_allele_elements[0]
                p = alternate_allele_elements[1]
            elif re.search("^.*\].*\]$", self.alternate_allele):                # t]p] reverse comp piece extending left of p is joined after t
                orientation = TranslocationOrientations.ORIENTATION_2
                alternate_allele_elements = self.alternate_allele.split(']')
                t = alternate_allele_elements[0]
                p = alternate_allele_elements[1]
            elif re.search("^\].*\].*$", self.alternate_allele):                # ]p]t piece extending to the left of p is joined before t
                orientation = TranslocationOrientations.ORIENTATION_3
                alternate_allele_elements = self.alternate_allele.split(']')
                t = alternate_allele_elements[2]
                p = alternate_allele_elements[1]
            elif re.search("^\[.*\[.*$", self.alternate_allele):                # [p[t  reverse comp piece extending right of p is joined before t
                orientation = TranslocationOrientations.ORIENTATION_4
                alternate_allele_elements = self.alternate_allele.split('[')
                t = alternate_allele_elements[2]
                p = alternate_allele_elements[1]
            else:
                raise Exception('Unknown ALT format to infer translocation orientation type: %s' % self.alternate_allele)

            if p == '%s:%i' % (self.chromosome_1, self.position_1):
                p_chromosome = self.chromosome_1
                p_position = self.position_1
                t_chromosome = self.chromosome_2
                t_position = self.position_2
            elif p == '%s:%i' % (self.chromosome_2, self.position_2):
                p_chromosome = self.chromosome_2
                p_position = self.position_2
                t_chromosome = self.chromosome_1
                t_position = self.position_1
            else:
                raise Exception('Positions for p and t could not be inferred from self.alternate_allele: %s' % self.alternate_allele)
            return orientation, t_chromosome, t_position, p_chromosome, p_position
        else:
            raise Exception('This VariantCall object does not encode a translocation. '
                            'Therefore translocation orientation cannot be inferred.')

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self.to_dataframe_row())

    def to_dataframe_row(self) -> Dict:
        data = {
            'variant_call_id': [self.id],
            'sample_id': [self.sample_id],
            'chromosome_1': [self.chromosome_1],
            'position_1': [self.position_1],
            'chromosome_2': [self.chromosome_2],
            'position_2': [self.position_2],
            'variant_type': [self.variant_type],
            'reference_allele': [self.reference_allele],
            'alternate_allele': [self.alternate_allele],
            'source_id': [self.source_id],
            'phase_block_id': [self.phase_block_id],
            'clone_id': [self.clone_id],
            'nucleic_acid': [self.nucleic_acid],
            'variant_calling_method': [self.variant_calling_method],
            'sequencing_platform': [self.sequencing_platform],
            'filter': [self.filter],
            'quality_score': [self.quality_score],
            'precise': [self.precise],
            'variant_subtype': [self.variant_subtype],
            'variant_size': [self.variant_size],
            'reference_allele_read_count': [self.reference_allele_read_count],
            'alternate_allele_read_count': [self.alternate_allele_read_count],
            'total_read_count': [self.total_read_count],
            'alternate_allele_fraction': [self.alternate_allele_fraction],
            'alternate_allele_read_ids': [';'.join(self.alternate_allele_read_ids)],
            'variant_sequences': [';'.join(self.variant_sequences)],
            'tags': [';'.join([str(i) for i in list(self.tags)])],
            'average_alignment_score_window': [self.average_alignment_score_window],
            'position_1_average_alignment_score': [self.position_1_average_alignment_score],
            'position_2_average_alignment_score': [self.position_2_average_alignment_score]
        }
        attributes = []
        for key, val in self.attributes.items():
            attributes.append('%s=%s' % (key, val))
        data['attributes'] = [';'.join(attributes)]

        pos_1_annotation_annotator = []
        pos_1_annotation_annotator_version = []
        pos_1_annotation_gene_id = []
        pos_1_annotation_gene_id_stable = []
        pos_1_annotation_gene_name = []
        pos_1_annotation_gene_strand = []
        pos_1_annotation_gene_type = []
        pos_1_annotation_gene_version = []
        pos_1_annotation_region = []
        pos_1_annotation_species = []
        for i in self.position_1_annotations:
            pos_1_annotation_annotator.append(i.annotator)
            pos_1_annotation_annotator_version.append(i.annotator_version)
            pos_1_annotation_gene_id.append(i.gene_id)
            pos_1_annotation_gene_id_stable.append(i.gene_id_stable)
            pos_1_annotation_gene_name.append(i.gene_name)
            pos_1_annotation_gene_strand.append(i.gene_strand)
            pos_1_annotation_gene_type.append(i.gene_type)
            pos_1_annotation_gene_version.append(i.gene_version)
            pos_1_annotation_region.append(i.region)
            pos_1_annotation_species.append(i.species)

        pos_2_annotation_annotator = []
        pos_2_annotation_annotator_version = []
        pos_2_annotation_gene_id = []
        pos_2_annotation_gene_id_stable = []
        pos_2_annotation_gene_name = []
        pos_2_annotation_gene_strand = []
        pos_2_annotation_gene_type = []
        pos_2_annotation_gene_version = []
        pos_2_annotation_region = []
        pos_2_annotation_species = []
        for i in self.position_2_annotations:
            pos_2_annotation_annotator.append(i.annotator)
            pos_2_annotation_annotator_version.append(i.annotator_version)
            pos_2_annotation_gene_id.append(i.gene_id)
            pos_2_annotation_gene_id_stable.append(i.gene_id_stable)
            pos_2_annotation_gene_name.append(i.gene_name)
            pos_2_annotation_gene_strand.append(i.gene_strand)
            pos_2_annotation_gene_type.append(i.gene_type)
            pos_2_annotation_gene_version.append(i.gene_version)
            pos_2_annotation_region.append(i.region)
            pos_2_annotation_species.append(i.species)

        data['position_1_annotation_annotator'] = [';'.join(pos_1_annotation_annotator)]
        data['position_1_annotation_annotator_version'] = [';'.join(pos_1_annotation_annotator_version)]
        data['position_1_annotation_gene_id'] = [';'.join(pos_1_annotation_gene_id)]
        data['position_1_annotation_gene_id_stable'] = [';'.join(pos_1_annotation_gene_id_stable)]
        data['position_1_annotation_gene_name'] = [';'.join(pos_1_annotation_gene_name)]
        data['position_1_annotation_gene_strand'] = [';'.join(pos_1_annotation_gene_strand)]
        data['position_1_annotation_gene_type'] = [';'.join(pos_1_annotation_gene_type)]
        data['position_1_annotation_gene_version'] = [';'.join(pos_1_annotation_gene_version)]
        data['position_1_annotation_region'] = [';'.join(pos_1_annotation_region)]
        data['position_1_annotation_species'] = [';'.join(pos_1_annotation_species)]

        data['position_2_annotation_annotator'] = [';'.join(pos_2_annotation_annotator)]
        data['position_2_annotation_annotator_version'] = [';'.join(pos_2_annotation_annotator_version)]
        data['position_2_annotation_gene_id'] = [';'.join(pos_2_annotation_gene_id)]
        data['position_2_annotation_gene_id_stable'] = [';'.join(pos_2_annotation_gene_id_stable)]
        data['position_2_annotation_gene_name'] = [';'.join(pos_2_annotation_gene_name)]
        data['position_2_annotation_gene_strand'] = [';'.join(pos_2_annotation_gene_strand)]
        data['position_2_annotation_gene_type'] = [';'.join(pos_2_annotation_gene_type)]
        data['position_2_annotation_gene_version'] = [';'.join(pos_2_annotation_gene_version)]
        data['position_2_annotation_region'] = [';'.join(pos_2_annotation_region)]
        data['position_2_annotation_species'] = [';'.join(pos_2_annotation_species)]

        return data

    def to_dict(self) -> Dict:
        data = {
            'id': self.id,
            'sample_id': self.sample_id,
            'chromosome_1': self.chromosome_1,
            'position_1': self.position_1,
            'chromosome_2': self.chromosome_2,
            'position_2': self.position_2,
            'variant_type': self.variant_type,
            'reference_allele': self.reference_allele,
            'alternate_allele': self.alternate_allele,
            'source_id': self.source_id,
            'phase_block_id': self.phase_block_id,
            'clone_id': self.clone_id,
            'nucleic_acid': self.nucleic_acid,
            'variant_calling_method': self.variant_calling_method,
            'sequencing_platform': self.sequencing_platform,
            'filter': self.filter,
            'quality_score': self.quality_score,
            'precise': self.precise,
            'variant_subtype': self.variant_subtype,
            'variant_size': self.variant_size,
            'total_read_count': self.total_read_count,
            'reference_allele_read_count': self.reference_allele_read_count,
            'alternate_allele_read_count': self.alternate_allele_read_count,
            'alternate_allele_fraction': self.alternate_allele_fraction,
            'alternate_allele_read_ids': list(self.alternate_allele_read_ids),
            'variant_sequences': list(self.variant_sequences),
            'tags': list(self.tags),
            'average_alignment_score_window': self.average_alignment_score_window,
            'position_1_average_alignment_score': self.position_1_average_alignment_score,
            'position_2_average_alignment_score': self.position_2_average_alignment_score
        }
        attributes = {}
        for key, value in self.attributes.items():
            attributes[key] = str(value)
        data['attributes'] = attributes
        data['position_1_annotations'] = [attribute.to_dict() for attribute in self.position_1_annotations]
        data['position_2_annotations'] = [attribute.to_dict() for attribute in self.position_2_annotations]
        return data
