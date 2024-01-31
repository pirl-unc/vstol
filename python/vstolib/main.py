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
The purpose of this python3 script is to implement main APIs.
"""


import copy
import pandas as pd
from collections import defaultdict
from typing import List, Tuple
from .annotator import Annotator
from .constants import VariantCallingMethods, VariantCallTags
from .default import *
from .genomic_ranges_list import GenomicRangesList
from .logging import get_logger
from .utilities import is_repeated_sequence
from .variant import Variant
from .variants_list import VariantsList
from .variant_filter import VariantFilter
from .vcf.cutesv import parse_cutesv_callset
from .vcf.dbsnp import parse_dbsnp_callset
from .vcf.deepvariant import parse_deepvariant_callset
from .vcf.delly2 import parse_delly2_somatic_callset
from .vcf.gatk4_mutect2 import parse_gatk4_mutect2_callset
from .vcf.lumpy import parse_lumpy_somatic_callset
from .vcf.pbsv import parse_pbsv_callset
from .vcf.sniffles2 import parse_sniffles2_callset
from .vcf.strelka2 import parse_strelka2_somatic_callset
from .vcf.svim import parse_svim_callset


logger = get_logger(__name__)


def annotate(
        variants_list: VariantsList,
        annotator: Annotator
) -> VariantsList:
    """
    Annotate a VariantsList object.

    Parameters:
        variants_list   :   VariantsList.
        annotator       :   Annotator.

    Returns:
        VariantsList
    """
    return annotator.annotate(variants_list=variants_list)


def diff(
        target_variants_list: VariantsList,
        query_variants_lists: List[VariantsList],
        max_neighbor_distance: int = DIFF_MAX_NEIGHBOR_DISTANCE,
        num_threads: int = NUM_THREADS,
        match_all_breakpoints: bool = DIFF_MATCH_ALL_BREAKPOINTS,
        match_variant_types: bool = DIFF_MATCH_VARIANT_TYPES
) -> VariantsList:
    """
    Diff variants lists.

    Parameters:
        target_variants_list    :   VariantsList object.
        query_variants_lists    :   List of VariantsList objects.
        max_neighbor_distance   :   Maximum neighbor distance.
        num_threads             :   Number of threads.
        match_all_breakpoints   :   If True, for two VariantCall objects
                                    to be considered an intersect, both
                                    pairs of breakpoints must match.
        match_variant_types     :   If True, for two VariantCall objects
                                    to be considered an intersect, the
                                    variant types must match.

    Returns:
        VariantsList (with VariantCalls private to target_variants_list)
    """
    variants_list_diff = copy.deepcopy(target_variants_list)
    for variants_list in query_variants_lists:
        variants_list_diff = variants_list_diff.diff(
            variants_list=variants_list,
            num_threads=num_threads,
            max_neighbor_distance=max_neighbor_distance,
            match_all_breakpoints=match_all_breakpoints,
            match_variant_types=match_variant_types
        )
    logger.info('%i variants and %i variant calls in the target VariantsList before diff.' %
                (len(target_variants_list.variants), len(target_variants_list.variant_call_ids)))
    logger.info('%i variants and %i variant calls in the target VariantsList after diff.' %
                (len(variants_list_diff.variants), len(variants_list_diff.variant_call_ids)))
    return variants_list_diff


def filter(
        variants_list: VariantsList,
        variant_filters: List[VariantFilter] = None,
        num_threads: int = NUM_THREADS
) -> Tuple[VariantsList, VariantsList]:
    """
    Filter out VariantCall objects.

    Parameters:
        variants_list               :   VariantsList object.
        variant_filters             :   List of VariantFilter objects.
        num_threads                 :   Number of threads.

    Returns:
        Tuple[variants_list_passed,variants_list_rejected]
    """
    logger.info('%i variants and %i variant calls in the original list before filtering.' %
                (len(variants_list.variants), len(variants_list.variant_call_ids)))

    rejected_variant_call_ids_dict = defaultdict(list)

    # Step 1. Filter out variants based on VariantFilter
    # key   = variant ID
    # value = reasons why the variant was rejected
    passed_variants_list = variants_list.filter(
        variant_filters=variant_filters,
        num_threads=num_threads,
    )
    passed_variants_ids = set(passed_variants_list.variant_ids)

    # Step 2. Create a passed VariantsList and a rejected VariantsList
    variants_list_passed = VariantsList()
    variants_list_rejected = VariantsList()
    for variant in variants_list.variants:
        if variant.id in passed_variants_ids:
            for i in range(0, variant.num_variant_calls):
                variant.variant_calls[i].tags.add(VariantCallTags.PASSED)
            variants_list_passed.add_variant(variant=variant)
        else:
            for i in range(0, variant.num_variant_calls):
                variant.variant_calls[i].tags.add(VariantCallTags.FAILED_FILTER)
            variants_list_rejected.add_variant(variant=variant)

    logger.info('%i variants and %i variant calls in the passed VariantsList after filtering.' %
                (len(variants_list_passed.variants), len(variants_list_passed.variant_call_ids)))
    logger.info('%i variants and %i variant calls in the rejected VariantsList after filtering.' %
                (len(variants_list_rejected.variants), len(variants_list_rejected.variant_call_ids)))

    return variants_list_passed, variants_list_rejected


def filter_excluded_regions(
        variants_list: VariantsList,
        excluded_regions_list: GenomicRangesList = None,
        excluded_regions_padding: int = FILTER_EXCLUDED_REGION_PADDING,
        num_threads: int = NUM_THREADS
) -> Tuple[VariantsList, VariantsList]:
    """
    Filter out VariantCalls in excluded regions.

    Parameters:
        variants_list               :   VariantsList object.
        excluded_regions_list       :   GenomicRangesList object of regions to exclude.
        excluded_regions_padding    :   Number of bases to pad each region to exclude.
        num_threads                 :   Number of threads.

    Returns:
        Tuple[variants_list_passed,variants_list_rejected]
    """
    logger.info('%i variants and %i variant calls in the original list before filtering.' %
                (len(variants_list.variants), len(variants_list.variant_call_ids)))

    # Step 1. Filter out VariantCall objects overlapping the excluded regions
    rejected_variant_call_ids = set()
    overlapping_variant_call_ids = variants_list.overlap(
        genomic_ranges_list=excluded_regions_list,
        padding=excluded_regions_padding,
        num_threads=num_threads
    )
    for variant_call_id, genomic_ranges in overlapping_variant_call_ids:
        rejected_variant_call_ids.add(variant_call_id)
    logger.info('%i variant calls are near excluded regions.' % len(overlapping_variant_call_ids))

    # Step 2. Create a passed VariantsList and a rejected VariantsList
    variants_list_passed = VariantsList()
    variants_list_rejected = VariantsList()
    for variant in variants_list.variants:
        variant_calls_passed = []
        variant_calls_rejected = []
        for variant_call in variant.variant_calls:
            if variant_call.id in rejected_variant_call_ids:
                variant_call.tags.add(VariantCallTags.NEARBY_EXCLUDED_REGION)
                variant_calls_rejected.append(variant_call)
            else:
                variant_call.tags.add(VariantCallTags.PASSED)
                variant_calls_passed.append(variant_call)

        variant_passed = Variant(id=variant.id)
        variant_rejected = Variant(id=variant.id)
        for variant_call in variant_calls_passed:
            variant_passed.add_variant_call(variant_call=variant_call)
        for variant_call in variant_calls_rejected:
            variant_rejected.add_variant_call(variant_call=variant_call)

        if variant_passed.num_variant_calls > 0:
            variants_list_passed.add_variant(variant=variant_passed)
        if variant_rejected.num_variant_calls > 0:
            variants_list_rejected.add_variant(variant=variant_rejected)

    logger.info('%i variants and %i variant calls in the passed VariantsList after filtering.' %
                (len(variants_list_passed.variants), len(variants_list_passed.variant_call_ids)))
    logger.info('%i variants and %i variant calls in the rejected VariantsList after filtering.' %
                (len(variants_list_rejected.variants), len(variants_list_rejected.variant_call_ids)))

    return variants_list_passed, variants_list_rejected


def filter_homopolymeric_variants(
        variants_list: VariantsList,
        reference_genome_fasta_file: str,
        homopolymer_length: int = FILTER_HOMOPOLYMER_LENGTH,
        num_threads: int = NUM_THREADS
) -> Tuple[VariantsList, VariantsList]:
    """
    Filter out homopolymeric VariantCall objects.

    Parameters:
        variants_list                   :   VariantsList object.
        reference_genome_fasta_file     :   Reference genome FASTA file.
        homopolymer_length              :   Homopolymer length.
        num_threads                     :   Number of threads.

    Returns:
        Tuple[variants_list_passed,variants_list_rejected]
    """
    logger.info('%i variants and %i variant calls in the original list before filtering.' %
                (len(variants_list.variants), len(variants_list.variant_call_ids)))

    # Step 1. Filter out variants in homopolymer regions
    flanking_sequences = variants_list.find_breakpoint_flanking_sequences(
        reference_genome_fasta_file=reference_genome_fasta_file,
        num_threads=num_threads,
        length=homopolymer_length
    )
    rejected_variant_call_ids = set()
    for variant_id, variant_call_id, position, left_seq, right_seq in flanking_sequences:
        if len(left_seq) == homopolymer_length and is_repeated_sequence(sequence=left_seq):
            rejected_variant_call_ids.add(variant_call_id)
        if len(right_seq) == homopolymer_length and is_repeated_sequence(sequence=right_seq):
            rejected_variant_call_ids.add(variant_call_id)
    logger.info('%i variant calls are homopolymeric variant calls.' % len(rejected_variant_call_ids))

    # Step 2. Create a passed VariantsList and a rejected VariantsList
    variants_list_passed = VariantsList()
    variants_list_rejected = VariantsList()
    for variant in variants_list.variants:
        variant_calls_passed = []
        variant_calls_rejected = []
        for variant_call in variant.variant_calls:
            if variant_call.id in rejected_variant_call_ids:
                variant_call.tags.add(VariantCallTags.HOMOPOLYMER_REGION)
                variant_calls_rejected.append(variant_call)
            else:
                variant_call.tags.add(VariantCallTags.PASSED)
                variant_calls_passed.append(variant_call)

        variant_passed = Variant(id=variant.id)
        variant_rejected = Variant(id=variant.id)
        for variant_call in variant_calls_passed:
            variant_passed.add_variant_call(variant_call=variant_call)
        for variant_call in variant_calls_rejected:
            variant_rejected.add_variant_call(variant_call=variant_call)

        if variant_passed.num_variant_calls > 0:
            variants_list_passed.add_variant(variant=variant_passed)
        if variant_rejected.num_variant_calls > 0:
            variants_list_rejected.add_variant(variant=variant_rejected)

    logger.info('%i variants and %i variant calls in the passed VariantsList after identifying homopolymeric variant calls.' %
                (len(variants_list_passed.variants), len(variants_list_passed.variant_call_ids)))
    logger.info('%i variants and %i variant calls in the rejected VariantsList after identifying homopolymeric variant calls.' %
                (len(variants_list_rejected.variants), len(variants_list_rejected.variant_call_ids)))

    return variants_list_passed, variants_list_rejected


def intersect(
        variants_lists: List[VariantsList],
        num_threads: int = NUM_THREADS,
        max_neighbor_distance: int = INTERSECT_MAX_NEIGHBOR_DISTANCE,
        match_all_breakpoints: bool = MERGE_MATCH_ALL_BREAKPOINTS,
        match_variant_types: bool = MERGE_MATCH_VARIANT_TYPES
) -> VariantsList:
    """
    Return intersecting VariantsList.

    Parameters:
        variants_lists          :   List of VariantsList objects.
        num_threads             :   Number of threads.
        max_neighbor_distance   :   Maximum neighbor distance.
        match_all_breakpoints       :   If True, for two VariantCall objects to be considered
                                        intersecting, all breakpoints must match or be near each other.
        match_variant_types         :   If True, for two VariantCall objects to be considered
                                        intersecting, their variant types must match.

    Returns:
        VariantsList
    """
    return VariantsList.intersect(
        variants_lists=variants_lists,
        num_threads=num_threads,
        max_neighbor_distance=max_neighbor_distance,
        match_all_breakpoints=match_all_breakpoints,
        match_variant_types=match_variant_types
    )


def merge(
        variants_lists: List[VariantsList],
        num_threads: int = NUM_THREADS,
        max_neighbor_distance: int = MERGE_MAX_NEIGHBOR_DISTANCE,
        match_all_breakpoints: bool = MERGE_MATCH_ALL_BREAKPOINTS,
        match_variant_types: bool = MERGE_MATCH_VARIANT_TYPES
) -> VariantsList:
    """
    Merge VariantsList objects into one.

    Parameters:
        variants_lists              :   List of VariantsList objects.
        num_threads                 :   Number of threads.
        max_neighbor_distance       :   Maximum neighbor distance.
        match_all_breakpoints       :   If True, for two VariantCall objects to be considered
                                        intersecting, all breakpoints must match or be near each other.
        match_variant_types         :   If True, for two VariantCall objects to be considered
                                        intersecting, their variant types must match.

    Returns:
        VariantsList
    """
    return VariantsList.merge(
        variants_lists=variants_lists,
        num_threads=num_threads,
        max_neighbor_distance=max_neighbor_distance,
        match_all_breakpoints=match_all_breakpoints,
        match_variant_types=match_variant_types
    )


def overlap(
        variants_list: VariantsList,
        genomic_ranges_list: GenomicRangesList,
        padding: int = OVERLAP_PADDING,
        num_threads: int = NUM_THREADS
) -> VariantsList:
    """
    Overlap variants list.

    Parameters:
        variants_list:          :   VariantsList object.
        genomic_ranges_list:    :   GenomicRangesList object.
        padding:                :   Padding (applied to each breakpoint for every variant call).
        num_threads:            :   Number of threads.

    Returns:
        VariantsList
    """
    # Step 1. Identify overlapping variants
    overlapping_variant_call_ids = variants_list.overlap(
        genomic_ranges_list=genomic_ranges_list,
        num_threads=num_threads,
        padding=padding
    )
    overlapping_variant_call_ids = [i[0] for i in overlapping_variant_call_ids]

    # Step 2. Prepare variants list to output
    variants_list_overlapping = VariantsList()
    for variant in variants_list.variants:
        variant_ = Variant(id=variant.id)
        for variant_call in variant.variant_calls:
            if variant_call.id in overlapping_variant_call_ids:
                variant_.add_variant_call(variant_call=variant_call)
        if variant_.num_variant_calls > 0:
            variants_list_overlapping.add_variant(variant=variant_)

    logger.info('%i variants and %i variant calls overlap' %
                (len(variants_list_overlapping.variant_ids),
                 len(variants_list_overlapping.variant_call_ids)))
    return variants_list_overlapping


def vcf2tsv(
        df_vcf: pd.DataFrame,
        source_id: str,
        variant_calling_method: str,
        sequencing_platform: str,
        case_id: str = '',
        control_id: str = ''
) -> VariantsList:
    """
    Convert a VCF Pandas DataFrame to a VariantsList object.

    Parameters:
        df_vcf                  :   Pandas DataFrame (read from vcf.common.read_vcf_file).
        source_id               :   Source ID (e.g. patient ID or cell line sample ID).
        variant_calling_method  :   Variant calling method.
        sequencing_platform     :   Sequencing platform.
        case_id                 :   Case ID (only necessary if variant_calling_method is 'strelka2-somatic').
        control_id              :   Control ID (only necessary if variant_calling_method is 'strelka2-somatic').

    Returns:
        VariantsList
    """
    if variant_calling_method == VariantCallingMethods.CUTESV:
        variants_list = parse_cutesv_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    elif variant_calling_method == VariantCallingMethods.DEEPVARIANT:
        variants_list = parse_deepvariant_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    elif variant_calling_method == VariantCallingMethods.DELLY2_SOMATIC:
        variants_list = parse_delly2_somatic_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    elif variant_calling_method == VariantCallingMethods.GATK4_MUTECT2:
        variants_list = parse_gatk4_mutect2_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    elif variant_calling_method == VariantCallingMethods.LUMPY_SOMATIC:
        variants_list = parse_lumpy_somatic_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    elif variant_calling_method == VariantCallingMethods.PBSV:
        variants_list = parse_pbsv_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    elif variant_calling_method == VariantCallingMethods.SNIFFLES2:
        variants_list = parse_sniffles2_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    elif variant_calling_method == VariantCallingMethods.STRELKA2_SOMATIC:
        variants_list = parse_strelka2_somatic_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id,
            case_id=case_id,
            control_id=control_id
        )
    elif variant_calling_method == VariantCallingMethods.SVIM:
        variants_list = parse_svim_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    elif variant_calling_method == VariantCallingMethods.DBSNP:
        variants_list = parse_dbsnp_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    else:
        raise Exception('Unsupported variant calling method: %s' % variant_calling_method)
    return variants_list

