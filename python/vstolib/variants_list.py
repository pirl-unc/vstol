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
The purpose of this python3 script is to implement the VariantsList dataclass.
"""


import json
import pandas as pd
import pysam
import multiprocessing as mp
from collections import defaultdict, OrderedDict
from dataclasses import dataclass, field
from functools import partial
from typing import Dict, List, Tuple, Type
from vstolib import vstolibrs
from .genomic_range import GenomicRange
from .genomic_ranges_list import GenomicRangesList
from .logging import get_logger
from .utilities import get_typed_value, is_gzipped, retrieve_from_dict
from .variant import Variant
from .variant_call_annotation import VariantCallAnnotation
from .variant_call import VariantCall
from .variant_filter import VariantFilter
from .vcf.common import get_attribute_types


logger = get_logger(__name__)


@dataclass
class VariantsList:
    variants: List[Variant] = field(default_factory=list)

    # Dictionary for bookkeeping purposes
    # key   = variant ID
    # value = index
    _variants_dict: Dict[str, int] = field(default_factory=dict)

    def __post_init__(self):
        for i in range(0, len(self.variants)):
            self._variants_dict[self.variants[i].id] = i

    @property
    def size(self) -> int:
        return len(self.variants)

    @property
    def variant_call_ids(self) -> List[str]:
        variant_call_ids = []
        for variant in self.variants:
            for variant_call in variant.variant_calls:
                variant_call_ids.append(variant_call.id)
        return variant_call_ids

    @property
    def variant_ids(self) -> List[str]:
        variant_ids = []
        for variant in self.variants:
            variant_ids.append(variant.id)
        return variant_ids

    def add_variant(self, variant: Variant):
        """
        Add a Variant object.
        """
        curr_len = len(self.variants)
        self._variants_dict[variant.id] = curr_len
        self.variants.append(variant)

    def filter(
            self,
            variant_filters: List[VariantFilter],
            num_threads: int
    ) -> 'VariantsList':
        """
        Filter variants by a list of VariantFilter objects and return a
        VariantsList with variants that pass all VariantFilter objects.
        """
        # Step 1. Serialize VariantsList object
        variants_list_serialized = json.dumps(self.to_dict())

        # Step 2. Serialize VariantFilter objects
        variant_filters_serialized = []
        for variant_filter in variant_filters:
            variant_filters_serialized.append(json.dumps(variant_filter.to_dict()))

        # Step 3. Filter VariantsList object
        json_str = vstolibrs.filter_variants_list(
            variants_list_serialized,
            variant_filters_serialized,
            num_threads
        )

        # Step 4. Deserialize filtered VariantsList object
        return VariantsList.load_serialized_json(json_str=json_str)

    def find_breakpoint_flanking_sequences(
            self,
            reference_genome_fasta_file: str,
            num_threads: int,
            length: int,
    ) -> List[Tuple[str,str,int,str,str]]:
        """
        Find flanking sequences of every breakpoint.

        Parameters:
            reference_genome_fasta_file :   Reference genome FASTA file.
            num_threads                 :   Number of threads.
            length                      :   Flanking sequence length.

        Returns:
            List[Tuple[variant_id,variant_call_id,position,left_sequence,right_sequence]]
        """
        pool = mp.Pool(processes=num_threads)
        func = partial(VariantsList.find_breakpoint_flanking_sequences_worker, reference_genome_fasta_file, length)
        flanking_sequences_lists = pool.map(func, self.variants)
        pool.close()
        flanking_sequences = []
        for flanking_sequences_list in flanking_sequences_lists:
            flanking_sequences.extend(flanking_sequences_list)
        return flanking_sequences

    def get_variant(self, variant_id: str) -> Variant:
        return self.variants[self._variants_dict[variant_id]]

    def overlap(
            self,
            genomic_ranges_list: GenomicRangesList,
            num_threads: int
    ) -> List[Tuple[str, List[GenomicRange]]]:
        """
        Identify overlapping variants based on a GenomicRangesList object.

        Parameters:
            genomic_ranges_list         :   GenomicRangesList object.
            num_threads                 :   Number of threads.

        Returns:
            List[Tuple[variant_call_id,List[GenomicRange]]]
        """
        # Step 1. Serialize VariantsList object
        variants_list_serialized = json.dumps(self.to_dict())

        # Step 2. Serialize GenomicRangesList object
        genomic_ranges_list_serialized = json.dumps(genomic_ranges_list.to_dict())

        # Step 3. Identify Variant objects that overlap GenomicRange objects
        nearby_variant_call_ids = vstolibrs.overlap_variant_calls(
            variants_list_serialized,
            genomic_ranges_list_serialized,
            num_threads
        )

        # Step 4. Get Variant and GenomicRange objects
        overlapping_variant_call_ids = []
        for variant_call_id, genomic_range_ids in nearby_variant_call_ids.items():
            genomic_ranges = []
            for genomic_range_id in genomic_range_ids:
                genomic_range = genomic_ranges_list.get_genomic_range(genomic_range_id)
                genomic_ranges.append(genomic_range)
            overlapping_variant_call_ids.append((variant_call_id, genomic_ranges))
        return overlapping_variant_call_ids

    def subtract(
            self,
            variants_list: 'VariantsList',
            num_threads: int,
            max_neighbor_distance: int,
            match_all_breakpoints: bool,
            match_variant_types: bool,
            min_ins_size_overlap: float,
            min_del_size_overlap: float
    ) -> 'VariantsList':
        """
        Subtract another variants list.

        Parameters:
            variants_list               :   VariantsList object.
            num_threads                 :   Number of threads.
            max_neighbor_distance:          Maximum neighbor distance.
            match_all_breakpoints       :   If True, for two VariantCall objects to be considered
                                            intersecting, all breakpoints must match or be near each other.
            match_variant_types         :   If True, for two VariantCall objects to be considered
                                            intersecting, their variant types must match.
            min_ins_size_overlap        :   Minimum insertion size overlap.
            min_del_size_overlap        :   Minimum deletion size overlap.

        Returns:
            VariantsList
        """
        # Step 1. Serialize VariantsList objects
        vl_a = json.dumps(self.to_dict())
        vl_b = json.dumps(variants_list.to_dict())

        # Step 2. Subtract B from A.
        json_str = vstolibrs.subtract_variants_list(
            vl_a,
            vl_b,
            num_threads,
            max_neighbor_distance,
            match_all_breakpoints,
            match_variant_types,
            min_ins_size_overlap,
            min_del_size_overlap
        )

        return VariantsList.load_serialized_json(json_str=json_str)

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self.to_dataframe_row())

    def to_dataframe_row(self) -> Dict:
        data = {
            'variant_id': [],
            'variant_call_id': [],
            'sample_id': [],
            'chromosome_1': [],
            'position_1': [],
            'chromosome_2': [],
            'position_2': [],
            'variant_type': [],
            'reference_allele': [],
            'alternate_allele': [],
            'source_id': [],
            'phase_block_id': [],
            'clone_id': [],
            'nucleic_acid': [],
            'variant_calling_method': [],
            'sequencing_platform': [],
            'filter': [],
            'quality_score': [],
            'precise': [],
            'variant_subtype': [],
            'variant_size': [],
            'reference_allele_read_count': [],
            'alternate_allele_read_count': [],
            'total_read_count': [],
            'alternate_allele_fraction': [],
            'alternate_allele_read_ids': [],
            'variant_sequences': [],
            'tags': [],
            'attributes': [],
            'average_alignment_score_window': [],
            'position_1_average_alignment_score': [],
            'position_2_average_alignment_score': [],
            'position_1_annotation_annotator': [],
            'position_1_annotation_annotator_version': [],
            'position_1_annotation_gene_id': [],
            'position_1_annotation_gene_id_stable': [],
            'position_1_annotation_gene_name': [],
            'position_1_annotation_gene_strand': [],
            'position_1_annotation_gene_type': [],
            'position_1_annotation_gene_version': [],
            'position_1_annotation_transcript_id': [],
            'position_1_annotation_transcript_id_stable': [],
            'position_1_annotation_transcript_name': [],
            'position_1_annotation_transcript_strand': [],
            'position_1_annotation_transcript_type': [],
            'position_1_annotation_transcript_version': [],
            'position_1_annotation_exon_id': [],
            'position_1_annotation_exon_id_stable': [],
            'position_1_annotation_region': [],
            'position_1_annotation_species': [],
            'position_2_annotation_annotator': [],
            'position_2_annotation_annotator_version': [],
            'position_2_annotation_gene_id': [],
            'position_2_annotation_gene_id_stable': [],
            'position_2_annotation_gene_name': [],
            'position_2_annotation_gene_strand': [],
            'position_2_annotation_gene_type': [],
            'position_2_annotation_gene_version': [],
            'position_2_annotation_transcript_id': [],
            'position_2_annotation_transcript_id_stable': [],
            'position_2_annotation_transcript_name': [],
            'position_2_annotation_transcript_strand': [],
            'position_2_annotation_transcript_type': [],
            'position_2_annotation_transcript_version': [],
            'position_2_annotation_exon_id': [],
            'position_2_annotation_exon_id_stable': [],
            'position_2_annotation_region': [],
            'position_2_annotation_species': []
        }
        for variant in self.variants:
            for key, values in variant.to_dataframe_row().items():
                for value in values:
                    data[key].append(value)
        return data

    def to_dict(self) -> Dict:
        data = {
            'variants': [variant.to_dict() for variant in self.variants]
        }
        return data

    @staticmethod
    def compare(
            a: 'VariantsList',
            b: 'VariantsList',
            num_threads: int,
            max_neighbor_distance: int,
            match_all_breakpoints: bool,
            match_variant_types: bool,
            min_ins_size_overlap: float,
            min_del_size_overlap: float
    ) -> Tuple['VariantsList','VariantsList','VariantsList']:
        """
        Compare two VariantsList objects and return three VariantsList objects: shared, specific to a, and specific to b.

        Parameters:
            a                       :   VariantsList object.
            b                       :   VariantsList object.
            num_threads             :   Number of threads.
            max_neighbor_distance   :   Maximum neighbor distance.
                                        This value is used to decide if
                                        a VariantCall should be appended to an existing Variant.
                                        If there exists a VariantCall in a given Variant
                                        where the distances to both position_1 and position_2
                                        are equal to or less than max_neighbor_distance
                                        to position_1 and position_2 of the specified VariantCall,
                                        respectively, then the specified VariantCall is
                                        appended to the Variant. If such VariantCall
                                        is not identified, then a new Variant is constructed
                                        and added to self.variants.
            match_all_breakpoints   :   If True, for two VariantCall objects to be considered
                                        intersecting, all breakpoints must match or be near each other.
            match_variant_types     :   If True, for two VariantCall objects to be considered
                                        intersecting, their variant types must match.
            min_ins_size_overlap    :   Minimum insertion size overlap.
            min_del_size_overlap    :   Minimum deletion size overlap.

        Returns:
            Tuples[shared VariantList, a-specific VariantList, b-specific VariantList]
        """
        # Step 1. Serialize all VariantsList objects
        variants_lists_serialized = []
        variants_lists_serialized.append(json.dumps(a.to_dict()))
        variants_lists_serialized.append(json.dumps(b.to_dict()))

        # Step 2. Compare VariantsList objects
        json_str = vstolibrs.compare_variants_lists(
            variants_lists_serialized,
            num_threads,
            max_neighbor_distance,
            match_all_breakpoints,
            match_variant_types,
            min_ins_size_overlap,
            min_del_size_overlap
        )

        variants_lists = VariantsList.load_serialized_json_list(json_str=json_str)
        assert(len(variants_lists) == 3)
        return variants_lists[0], variants_lists[1], variants_lists[2]

    @staticmethod
    def find_breakpoint_flanking_sequences_worker(
            reference_genome_fasta_file: str,
            length: int,
            variant: Variant
    ) -> List[Tuple[str,str,int,str,str]]:
        fasta = pysam.FastaFile(reference_genome_fasta_file)
        flanking_sequences = []
        for variant_call in variant.variant_calls:
            left_1 = right_1 = left_2 = right_2 = ''
            try:
                left_1 = fasta.fetch(variant_call.chromosome_1, variant_call.position_1 - length, variant_call.position_1)
                right_1 = fasta.fetch(variant_call.chromosome_1, variant_call.position_1 - 1, variant_call.position_1 + length - 1)
                left_2 = fasta.fetch(variant_call.chromosome_2, variant_call.position_2 - length, variant_call.position_2)
                right_2 = fasta.fetch(variant_call.chromosome_2, variant_call.position_2 - 1, variant_call.position_2 + length - 1)
            except:
                pass
            flanking_sequences.append(
                (variant.id, variant_call.id, variant_call.position_1, left_1, right_1)
            )
            flanking_sequences.append(
                (variant.id, variant_call.id, variant_call.position_2, left_2, right_2)
            )
        return flanking_sequences

    @staticmethod
    def intersect(
            variants_lists: List['VariantsList'],
            num_threads: int,
            max_neighbor_distance: int,
            match_all_breakpoints: bool,
            match_variant_types: bool,
            min_ins_size_overlap: float,
            min_del_size_overlap: float
    ) -> 'VariantsList':
        """
        Return an intersecting VariantsList.

        Parameters:
            variants_lists                  :   List of VariantsList objects.
            num_threads                     :   Number of threads.
            max_neighbor_distance           :   Maximum neighbor distance.
                                                This value is used to decide if
                                                a VariantCall should be appended to an existing Variant.
                                                If there exists a VariantCall in a given Variant
                                                where the distances to both position_1 and position_2
                                                are equal to or less than max_neighbor_distance
                                                to position_1 and position_2 of the specified VariantCall,
                                                respectively, then the specified VariantCall is
                                                appended to the Variant. If such VariantCall
                                                is not identified, then a new Variant is constructed
                                                and added to self.variants.
            match_all_breakpoints           :   If True, for two VariantCall objects to be considered
                                                intersecting, all breakpoints must match or be near each other.
            match_variant_types             :   If True, for two VariantCall objects to be considered
                                                intersecting, their variant types must match.
            min_ins_size_overlap            :   Minimum insertion size overlap.
            min_del_size_overlap            :   Minimum deletion size overlap.

        Returns:
            VariantList
        """
        # Step 1. Serialize all VariantsList objects
        variants_lists_serialized = []
        for variants_list in variants_lists:
            variants_lists_serialized.append(json.dumps(variants_list.to_dict()))

        # Step 2. Intersect VariantsList objects
        json_str = vstolibrs.intersect_variants_lists(
            variants_lists_serialized,
            num_threads,
            max_neighbor_distance,
            match_all_breakpoints,
            match_variant_types,
            min_ins_size_overlap,
            min_del_size_overlap
        )

        return VariantsList.load_serialized_json(json_str=json_str)

    @staticmethod
    def load_dataframe(df: pd.DataFrame) -> 'VariantsList':
        """
        Load a Pandas DataFrame and return a VariantsList object.

        Parameters:
            df      :   Pandas DataFrame. Mandatory headers:
                        'variant_id', 'variant_call_id', 'sample_id',
                        'chromosome_1', 'position_1', 'chromosome_2',
                        'position_2', 'variant_type', 'reference_allele',
                        'alternate_allele'.

        Returns:
            VariantsList
        """
        columns = df.columns.values.tolist()
        if 'variant_id' not in columns:
            raise Exception("The column 'variant_id' must exist.")
        if 'variant_call_id' not in columns:
            raise Exception("The column 'variant_call_id' must exist.")
        if 'sample_id' not in columns:
            raise Exception("The column 'sample_id' must exist.")
        if 'chromosome_1' not in columns:
            raise Exception("The column 'chromosome_1' must exist.")
        if 'position_1' not in columns:
            raise Exception("The column 'position_1' must exist.")
        if 'chromosome_2' not in columns:
            raise Exception("The column 'chromosome_2' must exist.")
        if 'position_2' not in columns:
            raise Exception("The column 'position_2' must exist.")
        if 'variant_type' not in columns:
            raise Exception("The column 'variant_type' must exist.")
        if 'reference_allele' not in columns:
            raise Exception("The column 'reference_allele' must exist.")
        if 'alternate_allele' not in columns:
            raise Exception("The column 'alternate_allele' must exist.")

        variants: Dict[str, Variant] = {}
        for row in df.to_dict('records'):
            # Mandatory fields
            variant_id = str(row['variant_id'])
            variant_call = VariantCall(
                id=str(row['variant_call_id']),
                sample_id=str(row['sample_id']),
                chromosome_1=str(row['chromosome_1']),
                position_1=int(row['position_1']),
                chromosome_2=str(row['chromosome_2']),
                position_2=int(row['position_2']),
                variant_type=str(row['variant_type']),
                reference_allele=str(row['reference_allele']),
                alternate_allele=str(row['alternate_allele'])
            )

            # Optional fields
            variant_call.source_id = retrieve_from_dict(dct=row, key='source_id', default_value='', type=str)
            variant_call.phase_block_id = retrieve_from_dict(dct=row, key='phase_block_id', default_value='', type=str)
            variant_call.clone_id = retrieve_from_dict(dct=row, key='clone_id', default_value='', type=str)
            variant_call.nucleic_acid = retrieve_from_dict(dct=row, key='nucleic_acid', default_value='', type=str)
            variant_call.variant_calling_method = retrieve_from_dict(dct=row, key='variant_calling_method', default_value='', type=str)
            variant_call.sequencing_platform = retrieve_from_dict(dct=row, key='sequencing_platform', default_value='', type=str)
            variant_call.filter = retrieve_from_dict(dct=row, key='filter', default_value='', type=str)
            variant_call.quality_score = retrieve_from_dict(dct=row, key='quality_score', default_value=-1.0, type=float)
            variant_call.precise = retrieve_from_dict(dct=row, key='precise', default_value='', type=str)
            variant_call.variant_subtype = retrieve_from_dict(dct=row, key='variant_subtype', default_value='', type=str)
            variant_call.variant_size = retrieve_from_dict(dct=row, key='variant_size', default_value=-1, type=int)
            variant_call.total_read_count = retrieve_from_dict(dct=row, key='total_read_count', default_value=-1, type=int)
            variant_call.reference_allele_read_count = retrieve_from_dict(dct=row, key='reference_allele_read_count', default_value=-1, type=int)
            variant_call.alternate_allele_read_count = retrieve_from_dict(dct=row, key='alternate_allele_read_count', default_value=-1, type=int)
            variant_call.alternate_allele_fraction = retrieve_from_dict(dct=row, key='alternate_allele_fraction', default_value=-1.0, type=float)
            variant_call.average_alignment_score_window = retrieve_from_dict(dct=row, key='average_alignment_score_window', default_value=-1, type=int)
            variant_call.position_1_average_alignment_score = retrieve_from_dict(dct=row, key='position_1_average_alignment_score', default_value=-1.0, type=float)
            variant_call.position_2_average_alignment_score = retrieve_from_dict(dct=row, key='position_2_average_alignment_score', default_value=-1.0, type=float)
            alternate_allele_read_ids = retrieve_from_dict(dct=row, key='alternate_allele_read_ids', default_value='', type=str)
            variant_sequences = retrieve_from_dict(dct=row, key='variant_sequences', default_value='', type=str)
            attributes = retrieve_from_dict(dct=row, key='attributes', default_value='', type=str)
            tags = retrieve_from_dict(dct=row, key='tags', default_value='', type=str)

            position_1_annotation_annotator = retrieve_from_dict(dct=row, key='position_1_annotation_annotator', default_value='', type=str)
            position_1_annotation_annotator_version = retrieve_from_dict(dct=row, key='position_1_annotation_annotator_version', default_value='', type=str)
            position_1_annotation_gene_id = retrieve_from_dict(dct=row, key='position_1_annotation_gene_id', default_value='', type=str)
            position_1_annotation_gene_id_stable = retrieve_from_dict(dct=row, key='position_1_annotation_gene_id_stable', default_value='', type=str)
            position_1_annotation_gene_name = retrieve_from_dict(dct=row, key='position_1_annotation_gene_name', default_value='', type=str)
            position_1_annotation_gene_strand = retrieve_from_dict(dct=row, key='position_1_annotation_gene_strand', default_value='', type=str)
            position_1_annotation_gene_type = retrieve_from_dict(dct=row, key='position_1_annotation_gene_type', default_value='', type=str)
            position_1_annotation_gene_version = retrieve_from_dict(dct=row, key='position_1_annotation_gene_version', default_value='', type=str)
            position_1_annotation_transcript_id = retrieve_from_dict(dct=row, key='position_1_annotation_transcript_id', default_value='', type=str)
            position_1_annotation_transcript_id_stable = retrieve_from_dict(dct=row, key='position_1_annotation_transcript_id_stable', default_value='', type=str)
            position_1_annotation_transcript_name = retrieve_from_dict(dct=row, key='position_1_annotation_transcript_name', default_value='', type=str)
            position_1_annotation_transcript_strand = retrieve_from_dict(dct=row, key='position_1_annotation_transcript_strand', default_value='', type=str)
            position_1_annotation_transcript_type = retrieve_from_dict(dct=row, key='position_1_annotation_transcript_type', default_value='', type=str)
            position_1_annotation_transcript_version = retrieve_from_dict(dct=row, key='position_1_annotation_transcript_version', default_value='', type=str)
            position_1_annotation_exon_id = retrieve_from_dict(dct=row, key='position_1_annotation_exon_id', default_value='', type=str)
            position_1_annotation_exon_id_stable = retrieve_from_dict(dct=row, key='position_1_annotation_exon_id_stable', default_value='', type=str)
            position_1_annotation_region = retrieve_from_dict(dct=row, key='position_1_annotation_region', default_value='', type=str)
            position_1_annotation_species = retrieve_from_dict(dct=row, key='position_1_annotation_species', default_value='', type=str)

            position_2_annotation_annotator = retrieve_from_dict(dct=row, key='position_2_annotation_annotator', default_value='', type=str)
            position_2_annotation_annotator_version = retrieve_from_dict(dct=row, key='position_2_annotation_annotator_version', default_value='', type=str)
            position_2_annotation_gene_id = retrieve_from_dict(dct=row, key='position_2_annotation_gene_id', default_value='', type=str)
            position_2_annotation_gene_id_stable = retrieve_from_dict(dct=row, key='position_2_annotation_gene_id_stable', default_value='', type=str)
            position_2_annotation_gene_name = retrieve_from_dict(dct=row, key='position_2_annotation_gene_name', default_value='', type=str)
            position_2_annotation_gene_strand = retrieve_from_dict(dct=row, key='position_2_annotation_gene_strand', default_value='', type=str)
            position_2_annotation_gene_type = retrieve_from_dict(dct=row, key='position_2_annotation_gene_type', default_value='', type=str)
            position_2_annotation_gene_version = retrieve_from_dict(dct=row, key='position_2_annotation_gene_version', default_value='', type=str)
            position_2_annotation_transcript_id = retrieve_from_dict(dct=row, key='position_2_annotation_transcript_id', default_value='', type=str)
            position_2_annotation_transcript_id_stable = retrieve_from_dict(dct=row, key='position_2_annotation_transcript_id_stable', default_value='', type=str)
            position_2_annotation_transcript_name = retrieve_from_dict(dct=row, key='position_2_annotation_transcript_name', default_value='', type=str)
            position_2_annotation_transcript_strand = retrieve_from_dict(dct=row, key='position_2_annotation_transcript_strand', default_value='', type=str)
            position_2_annotation_transcript_type = retrieve_from_dict(dct=row, key='position_2_annotation_transcript_type', default_value='', type=str)
            position_2_annotation_transcript_version = retrieve_from_dict(dct=row, key='position_2_annotation_transcript_version', default_value='', type=str)
            position_2_annotation_exon_id = retrieve_from_dict(dct=row, key='position_2_annotation_exon_id', default_value='', type=str)
            position_2_annotation_exon_id_stable = retrieve_from_dict(dct=row, key='position_2_annotation_exon_id_stable', default_value='', type=str)
            position_2_annotation_region = retrieve_from_dict(dct=row, key='position_2_annotation_region', default_value='', type=str)
            position_2_annotation_species = retrieve_from_dict(dct=row, key='position_2_annotation_species', default_value='', type=str)

            # Alternate allele read IDs
            if alternate_allele_read_ids != '':
                for read_id in alternate_allele_read_ids.split(';'):
                    variant_call.alternate_allele_read_ids.add(str(read_id))

            # Variant sequences
            if variant_sequences != '':
                for seq in variant_sequences.split(';'):
                    variant_call.variant_sequences.add(str(seq))

            # Attributes
            if attributes != '' and variant_call.variant_calling_method != '':
                attribute_types = get_attribute_types(variant_calling_method=variant_call.variant_calling_method)
                for attribute in attributes.split(';'):
                    attribute_key = attribute.split('=')[0]
                    attribute_value = attribute.split('=')[1]
                    if attribute_types[attribute_key] == int:
                        default_value = -1
                    elif attribute_types[attribute_key] == float:
                        default_value = -1.0
                    elif attribute_types[attribute_key] == str:
                        default_value = ''
                    elif attribute_types[attribute_key] == bool:
                        default_value = False
                    else:
                        raise Exception('Unknown variable type for %s' % attribute_types[attribute_key])
                    attribute_value = get_typed_value(
                        value=attribute_value,
                        default_value=default_value,
                        type=attribute_types[attribute_key]
                    )
                    if attribute_value != default_value:
                        variant_call.attributes[attribute_key] = attribute_value

            # Tags
            if tags != '':
                for tag in tags.split(';'):
                    variant_call.tags.add(str(tag))

            # Annotations
            if position_1_annotation_annotator != '':
                position_1_annotation_annotator = position_1_annotation_annotator.split(';')
                position_1_annotation_annotator_version = position_1_annotation_annotator_version.split(';')
                position_1_annotation_region = position_1_annotation_region.split(';')
                position_1_annotation_species = position_1_annotation_species.split(';')
                # Gene
                position_1_annotation_gene_id = position_1_annotation_gene_id.split(';')
                position_1_annotation_gene_id_stable = position_1_annotation_gene_id_stable.split(';')
                position_1_annotation_gene_name = position_1_annotation_gene_name.split(';')
                position_1_annotation_gene_strand = position_1_annotation_gene_strand.split(';')
                position_1_annotation_gene_type = position_1_annotation_gene_type.split(';')
                position_1_annotation_gene_version = position_1_annotation_gene_version.split(';')
                # Transcript
                position_1_annotation_transcript_id = position_1_annotation_transcript_id.split(';')
                position_1_annotation_transcript_id_stable = position_1_annotation_transcript_id_stable.split(';')
                position_1_annotation_transcript_name = position_1_annotation_transcript_name.split(';')
                position_1_annotation_transcript_strand = position_1_annotation_transcript_strand.split(';')
                position_1_annotation_transcript_type = position_1_annotation_transcript_type.split(';')
                position_1_annotation_transcript_version = position_1_annotation_transcript_version.split(';')
                # Exon
                position_1_annotation_exon_id = position_1_annotation_exon_id.split(';')
                position_1_annotation_exon_id_stable = position_1_annotation_exon_id_stable.split(';')
                for i in range(0, len(position_1_annotation_annotator)):
                    variant_call_annotation = VariantCallAnnotation(
                        annotator=position_1_annotation_annotator[i],
                        annotator_version=position_1_annotation_annotator_version[i],
                        region=position_1_annotation_region[i],
                        species=position_1_annotation_species[i],
                        gene_id=position_1_annotation_gene_id[i],
                        gene_id_stable=position_1_annotation_gene_id_stable[i],
                        gene_name=position_1_annotation_gene_name[i],
                        gene_strand=position_1_annotation_gene_strand[i],
                        gene_type=position_1_annotation_gene_type[i],
                        gene_version=position_1_annotation_gene_version[i],
                        transcript_id=position_1_annotation_transcript_id[i],
                        transcript_id_stable=position_1_annotation_transcript_id_stable[i],
                        transcript_name=position_1_annotation_transcript_name[i],
                        transcript_strand=position_1_annotation_transcript_strand[i],
                        transcript_type=position_1_annotation_transcript_type[i],
                        transcript_version=position_1_annotation_transcript_version[i],
                        exon_id=position_1_annotation_exon_id[i],
                        exon_id_stable=position_1_annotation_exon_id_stable[i]
                    )
                    variant_call.add_position_1_annotation(variant_call_annotation=variant_call_annotation)

            if position_2_annotation_annotator != '':
                position_2_annotation_annotator = position_2_annotation_annotator.split(';')
                position_2_annotation_annotator_version = position_2_annotation_annotator_version.split(';')
                position_2_annotation_region = position_2_annotation_region.split(';')
                position_2_annotation_species = position_2_annotation_species.split(';')
                # Gene
                position_2_annotation_gene_id = position_2_annotation_gene_id.split(';')
                position_2_annotation_gene_id_stable = position_2_annotation_gene_id_stable.split(';')
                position_2_annotation_gene_name = position_2_annotation_gene_name.split(';')
                position_2_annotation_gene_strand = position_2_annotation_gene_strand.split(';')
                position_2_annotation_gene_type = position_2_annotation_gene_type.split(';')
                position_2_annotation_gene_version = position_2_annotation_gene_version.split(';')
                # Transcript
                position_2_annotation_transcript_id = position_2_annotation_transcript_id.split(';')
                position_2_annotation_transcript_id_stable = position_2_annotation_transcript_id_stable.split(';')
                position_2_annotation_transcript_name = position_2_annotation_transcript_name.split(';')
                position_2_annotation_transcript_strand = position_2_annotation_transcript_strand.split(';')
                position_2_annotation_transcript_type = position_2_annotation_transcript_type.split(';')
                position_2_annotation_transcript_version = position_2_annotation_transcript_version.split(';')
                # Exon
                position_2_annotation_exon_id = position_2_annotation_exon_id.split(';')
                position_2_annotation_exon_id_stable = position_2_annotation_exon_id_stable.split(';')
                for i in range(0, len(position_2_annotation_annotator)):
                    variant_call_annotation = VariantCallAnnotation(
                        annotator=position_2_annotation_annotator[i],
                        annotator_version=position_2_annotation_annotator_version[i],
                        region=position_2_annotation_region[i],
                        species=position_2_annotation_species[i],
                        gene_id=position_2_annotation_gene_id[i],
                        gene_id_stable=position_2_annotation_gene_id_stable[i],
                        gene_name=position_2_annotation_gene_name[i],
                        gene_strand=position_2_annotation_gene_strand[i],
                        gene_type=position_2_annotation_gene_type[i],
                        gene_version=position_2_annotation_gene_version[i],
                        transcript_id=position_2_annotation_transcript_id[i],
                        transcript_id_stable=position_2_annotation_transcript_id_stable[i],
                        transcript_name=position_2_annotation_transcript_name[i],
                        transcript_strand=position_2_annotation_transcript_strand[i],
                        transcript_type=position_2_annotation_transcript_type[i],
                        transcript_version=position_2_annotation_transcript_version[i],
                        exon_id=position_2_annotation_exon_id[i],
                        exon_id_stable=position_2_annotation_exon_id_stable[i]
                    )
                    variant_call.add_position_2_annotation(variant_call_annotation=variant_call_annotation)

            if variant_id not in variants:
                variants[variant_id] = Variant(id=variant_id)
            variants[variant_id].add_variant_call(variant_call=variant_call)

        variants_list = VariantsList(variants=list(variants.values()))
        logger.info("Loaded %i variants and %i variant calls." % (variants_list.size, len(variants_list.variant_call_ids)))
        return variants_list

    @staticmethod
    def load_serialized_json(json_str: str) -> 'VariantsList':
        """
        Load a VariantsList object from a serialized JSON string.

        Parameters:
            json_str        :   JSON string.

        Returns:
            VariantsList
        """
        variants_list = VariantsList()
        variants_list_dict = json.loads(json_str)
        for variant_dict in variants_list_dict['variants']:
            variant = Variant(id=variant_dict['id'])
            for variant_call_dict in variant_dict['variant_calls']:
                position_1_annotations_dict = variant_call_dict['position_1_annotations']
                position_2_annotations_dict = variant_call_dict['position_2_annotations']
                del variant_call_dict['position_1_annotations']
                del variant_call_dict['position_2_annotations']
                variant_call = VariantCall(**variant_call_dict)
                for position_1_annotation_dict in position_1_annotations_dict:
                    variant_call.add_position_1_annotation(
                        variant_call_annotation=VariantCallAnnotation(**position_1_annotation_dict)
                    )
                for position_2_annotation_dict in position_2_annotations_dict:
                    variant_call.add_position_2_annotation(
                        variant_call_annotation=VariantCallAnnotation(**position_2_annotation_dict)
                    )
                variant.add_variant_call(variant_call=variant_call)
            variants_list.add_variant(variant=variant)
        return variants_list

    @staticmethod
    def load_serialized_json_list(json_str: str) -> List['VariantsList']:
        """
        Load a list of VariantsList objects from a serialized JSON string.

        Parameters:
            json_str        :   JSON string.

        Returns:
            List[VariantsList]
        """
        variants_lists = []
        for variants_list_dict in json.loads(json_str):
            variants_list = VariantsList()
            for variant_dict in variants_list_dict['variants']:
                variant = Variant(id=variant_dict['id'])
                for variant_call_dict in variant_dict['variant_calls']:
                    position_1_annotations_dict = variant_call_dict['position_1_annotations']
                    position_2_annotations_dict = variant_call_dict['position_2_annotations']
                    del variant_call_dict['position_1_annotations']
                    del variant_call_dict['position_2_annotations']
                    variant_call = VariantCall(**variant_call_dict)
                    for position_1_annotation_dict in position_1_annotations_dict:
                        variant_call.add_position_1_annotation(
                            variant_call_annotation=VariantCallAnnotation(**position_1_annotation_dict)
                        )
                    for position_2_annotation_dict in position_2_annotations_dict:
                        variant_call.add_position_2_annotation(
                            variant_call_annotation=VariantCallAnnotation(**position_2_annotation_dict)
                        )
                    variant.add_variant_call(variant_call=variant_call)
                variants_list.add_variant(variant=variant)
            variants_lists.append(variants_list)
        return variants_lists

    @staticmethod
    def merge(
            variants_lists: List['VariantsList'],
            num_threads: int,
            max_neighbor_distance: int,
            match_all_breakpoints: bool,
            match_variant_types: bool,
            min_ins_size_overlap: float,
            min_del_size_overlap: float
    ) -> 'VariantsList':
        """
        Merge a list of VariantsList objects and return a VariantsList object.

        Parameters:
            variants_lists                  :   List of VariantsList objects.
            num_threads                     :   Number of threads.
            max_neighbor_distance           :   Maximum neighbor distance.
                                                This value is used to decide if
                                                a VariantCall should be appended to an existing Variant.
                                                If there exists a VariantCall in a given Variant
                                                where the distances to both position_1 and position_2
                                                are equal to or less than max_neighbor_distance
                                                to position_1 and position_2 of the specified VariantCall,
                                                respectively, then the specified VariantCall is
                                                appended to the Variant. If such VariantCall
                                                is not identified, then a new Variant is constructed
                                                and added to self.variants.
            match_all_breakpoints           :   If True, for two VariantCall objects to be considered
                                                intersecting, all breakpoints must match or be near each other.
            match_variant_types             :   If True, for two VariantCall objects to be considered
                                                intersecting, their variant types must match.

        Returns:
            VariantList
        """
        # Step 1. Serialize all VariantsList objects
        variants_lists_serialized = []
        for variants_list in variants_lists:
            variants_lists_serialized.append(json.dumps(variants_list.to_dict()))

        # Step 2. Merge VariantsList objects
        json_str = vstolibrs.merge_variants_lists(
            variants_lists_serialized,
            num_threads,
            max_neighbor_distance,
            match_all_breakpoints,
            match_variant_types,
            min_ins_size_overlap,
            min_del_size_overlap
        )

        return VariantsList.load_serialized_json(json_str=json_str)

    @staticmethod
    def read_tsv_file(
            tsv_file: str,
            low_memory: bool = False,
            memory_map: bool = True
    ) -> 'VariantsList':
        """
        Read a TSV file and return a VariantsList object.

        Parameters:
            tsv_file    :   TSV file.
            low_memory  :   Low memory (default: False).
            memory_map  :   Map memory (default: True).

        Returns:
            VariantsList
        """
        if is_gzipped(tsv_file):
            df = pd.read_csv(tsv_file,
                             sep='\t',
                             compression='gzip',
                             low_memory=low_memory,
                             memory_map=memory_map)
        else:
            df = pd.read_csv(tsv_file,
                             sep='\t',
                             low_memory=low_memory,
                             memory_map=memory_map)
        return VariantsList.load_dataframe(df=df)

