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
import multiprocessing as mp
import numpy as np
import pandas as pd
import pysam
from collections import defaultdict
from functools import partial
from typing import List, Literal, Tuple
from vstolib import vstolibrs
from .annotator import Annotator
from .constants import CollapseStrategies, VariantCallingMethods
from .default import *
from .genomic_ranges_list import GenomicRangesList
from .logging import get_logger
from .utilities import is_repeated_sequence
from .variant import Variant
from .variants_list import VariantsList
from .variant_filter import VariantFilter
from .visualization import visualize
from .vcf.clairs import parse_clairs_callset
from .vcf.cutesv import parse_cutesv_callset
from .vcf.dbsnp import parse_dbsnp_callset
from .vcf.deepvariant import parse_deepvariant_callset
from .vcf.delly2 import parse_delly2_somatic_callset
from .vcf.gatk4_mutect2 import parse_gatk4_mutect2_callset
from .vcf.lumpy import parse_lumpy_somatic_callset
from .vcf.manta import parse_manta_somatic_callset
from .vcf.pbsv import parse_pbsv_callset
from .vcf.savana import parse_savana_callset
from .vcf.severus import parse_severus_callset
from .vcf.sniffles2 import parse_sniffles2_callset
from .vcf.strelka2 import parse_strelka2_somatic_callset
from .vcf.svim import parse_svim_callset
from .vcf.svisionpro import parse_svisionpro_callset


logger = get_logger(__name__)


def annotate(
        variants_list: VariantsList,
        annotator: Annotator,
        num_threads: int
) -> VariantsList:
    """
    Annotate a VariantsList object.

    Parameters:
        variants_list   :   VariantsList.
        annotator       :   Annotator.
        num_threads     :   Number of threads.

    Returns:
        VariantsList
    """
    return annotator.annotate(variants_list=variants_list,
                              num_processes=num_threads)


def collapse(
        variants_list: VariantsList,
        sample_id: str,
        strategy: str = Literal[CollapseStrategies.MAX_ALTERNATE_ALLELE_READ_COUNT]
) -> VariantsList:
    """
    Collapses (summarizes) a VariantsList such that each Variant has 1
    VariantCall.

    Args:
        variants_list   :   VariantsList object.
        sample_id       :   Sample ID to retain.
        strategy        :   Strategy (options: 'max_alternate_allele_read_count').

    Returns:
        VariantsList
    """
    if strategy == CollapseStrategies.MAX_ALTERNATE_ALLELE_READ_COUNT:
        variants_list_collapsed = VariantsList()
        for variant in variants_list.variants:
            target_idx = -1
            min_reads = -1
            matches_sample_id = False
            for i in range(0, len(variant.variant_calls)):
                variant_call = variant.variant_calls[i]
                if variant_call.sample_id == sample_id:
                    matches_sample_id = True
                    if variant_call.alternate_allele_read_count > min_reads:
                        target_idx = i
                        min_reads = variant_call.alternate_allele_read_count
            if matches_sample_id:
                if target_idx == -1:
                    target_idx = 0
                variant_ = Variant(id=variant.id)
                variant_.add_variant_call(variant_call=variant.variant_calls[target_idx])
                variants_list_collapsed.add_variant(variant=variant_)
    else:
        raise Exception('Unknown collapse strategy: %s' % strategy)
    logger.info("%i variants and %i variant calls in the collapsed VariantsList"
                % (variants_list_collapsed.size,
                   len(variants_list_collapsed.variant_call_ids)))
    return variants_list_collapsed


def compare(
        a: VariantsList,
        b: VariantsList,
        num_threads: int = NUM_THREADS,
        max_neighbor_distance: int = MAX_NEIGHBOR_DISTANCE,
        match_all_breakpoints: bool = MATCH_ALL_BREAKPOINTS,
        match_variant_types: bool = MATCH_VARIANT_TYPES,
        min_ins_size_overlap: float = MIN_INS_SIZE_OVERLAP,
        min_del_size_overlap: float = MIN_DEL_SIZE_OVERLAP
) -> Tuple[VariantsList, VariantsList, VariantsList]:
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
        min_ins_size_overlap        :   Minimum insertion size overlap.
        min_del_size_overlap        :   Minimum deletion size overlap.

    Returns:
        Tuple[shared VariantsList, a-specific VariantsList, b-specific VariantsList]
    """
    return VariantsList.compare(
        a=a,
        b=b,
        num_threads=num_threads,
        max_neighbor_distance=max_neighbor_distance,
        match_all_breakpoints=match_all_breakpoints,
        match_variant_types=match_variant_types,
        min_ins_size_overlap=min_ins_size_overlap,
        min_del_size_overlap=min_del_size_overlap
    )


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
            variants_list_passed.add_variant(variant=variant)
        else:
            variants_list_rejected.add_variant(variant=variant)

    logger.info('%i variants and %i variant calls in the passed VariantsList after filtering.' %
                (len(variants_list_passed.variants), len(variants_list_passed.variant_call_ids)))
    logger.info('%i variants and %i variant calls in the rejected VariantsList after filtering.' %
                (len(variants_list_rejected.variants), len(variants_list_rejected.variant_call_ids)))

    return variants_list_passed, variants_list_rejected


def filter_excluded_regions(
        variants_list: VariantsList,
        excluded_regions_list: GenomicRangesList,
        num_threads: int = NUM_THREADS
) -> Tuple[VariantsList, VariantsList]:
    """
    Filter out VariantCalls in excluded regions.

    Parameters:
        variants_list               :   VariantsList object.
        excluded_regions_list       :   GenomicRangesList object of regions to exclude.
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
                variant_calls_rejected.append(variant_call)
            else:
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
        homopolymer_length: int = HOMOPOLYMER_LENGTH,
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
                variant_calls_rejected.append(variant_call)
            else:
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
        max_neighbor_distance: int = MAX_NEIGHBOR_DISTANCE,
        match_all_breakpoints: bool = MATCH_ALL_BREAKPOINTS,
        match_variant_types: bool = MATCH_VARIANT_TYPES,
        min_ins_size_overlap: float = MIN_INS_SIZE_OVERLAP,
        min_del_size_overlap: float = MIN_DEL_SIZE_OVERLAP
) -> VariantsList:
    """
    Return intersecting VariantsList.

    Parameters:
        variants_lists              :   List of VariantsList objects.
        num_threads                 :   Number of threads.
        max_neighbor_distance       :   Maximum neighbor distance.
        match_all_breakpoints       :   If True, for two VariantCall objects to be considered
                                        intersecting, all breakpoints must match or be near each other.
        match_variant_types         :   If True, for two VariantCall objects to be considered
                                        intersecting, their variant types must match.
        min_ins_size_overlap        :   Minimum insertion size overlap.
        min_del_size_overlap        :   Minimum deletion size overlap.

    Returns:
        VariantsList
    """
    return VariantsList.intersect(
        variants_lists=variants_lists,
        num_threads=num_threads,
        max_neighbor_distance=max_neighbor_distance,
        match_all_breakpoints=match_all_breakpoints,
        match_variant_types=match_variant_types,
        min_ins_size_overlap=min_ins_size_overlap,
        min_del_size_overlap=min_del_size_overlap
    )


def merge(
        variants_lists: List[VariantsList],
        num_threads: int = NUM_THREADS,
        max_neighbor_distance: int = MAX_NEIGHBOR_DISTANCE,
        match_all_breakpoints: bool = MATCH_ALL_BREAKPOINTS,
        match_variant_types: bool = MATCH_VARIANT_TYPES,
        min_ins_size_overlap: float = MIN_INS_SIZE_OVERLAP,
        min_del_size_overlap: float = MIN_DEL_SIZE_OVERLAP
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
        min_ins_size_overlap        :   Minimum insertion size overlap.
        min_del_size_overlap        :   Minimum deletion size overlap.

    Returns:
        VariantsList
    """
    return VariantsList.merge(
        variants_lists=variants_lists,
        num_threads=num_threads,
        max_neighbor_distance=max_neighbor_distance,
        match_all_breakpoints=match_all_breakpoints,
        match_variant_types=match_variant_types,
        min_ins_size_overlap=min_ins_size_overlap,
        min_del_size_overlap=min_del_size_overlap
    )


def overlap(
        variants_list: VariantsList,
        genomic_ranges_list: GenomicRangesList,
        num_threads: int = NUM_THREADS
) -> VariantsList:
    """
    Overlap variants list.

    Parameters:
        variants_list:          :   VariantsList object.
        genomic_ranges_list:    :   GenomicRangesList object.
        num_threads:            :   Number of threads.

    Returns:
        VariantsList
    """
    # Step 1. Identify overlapping variants
    overlapping_variant_call_ids = variants_list.overlap(
        genomic_ranges_list=genomic_ranges_list,
        num_threads=num_threads
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


def score(
        variants_list: VariantsList,
        bam_file: str,
        window: int = WINDOW,
        num_threads: int = NUM_THREADS
) -> VariantsList:
    """
    Calculates average alignment score for each breakpoint.

    Args:
        variants_list   :   VariantsList object.
        bam_file        :   BAM file.
        window          :   Window (will be applied both upstream and downstream).
        num_threads     :   Number of threads.

    Returns:
        VariantsList
    """
    # Step 1. Get the regions
    regions = []
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    for variant in variants_list.variants:
        for variant_call in variant.variant_calls:
            # Position 1
            chromosome_length = bamfile.get_reference_length(variant_call.chromosome_1)
            start = variant_call.position_1 - window
            end = variant_call.position_1 + window
            if start < 0:
                start = 0
            if end > chromosome_length:
                end = chromosome_length
            regions.append((variant_call.chromosome_1,start,end))

            # Position 2
            chromosome_length = bamfile.get_reference_length(variant_call.chromosome_2)
            start = variant_call.position_2 - window
            end = variant_call.position_2 + window
            if start < 0:
                start = 0
            if end > chromosome_length:
                end = chromosome_length
            regions.append((variant_call.chromosome_2,start,end))

    # Step 2. Calculate the average alignment scores
    if len(regions) > 0:
        regions_scores = vstolibrs.calculate_average_alignment_scores(
            bam_file=bam_file,
            regions=regions,
            num_threads=num_threads
        )
        for variant in variants_list.variants:
            for variant_call in variant.variant_calls:
                # Position 1
                chromosome_length = bamfile.get_reference_length(variant_call.chromosome_1)
                start = variant_call.position_1 - window
                end = variant_call.position_1 + window
                if start < 0:
                    start = 0
                if end > chromosome_length:
                    end = chromosome_length
                variant_call.position_1_average_alignment_score = regions_scores[(variant_call.chromosome_1,start,end)]

                # Position 2
                chromosome_length = bamfile.get_reference_length(variant_call.chromosome_2)
                start = variant_call.position_2 - window
                end = variant_call.position_2 + window
                if start < 0:
                    start = 0
                if end > chromosome_length:
                    end = chromosome_length
                variant_call.position_2_average_alignment_score = regions_scores[(variant_call.chromosome_2,start,end)]

                variant_call.average_alignment_score_window = window

    return variants_list


def subtract(
        target_variants_list: VariantsList,
        query_variants_lists: List[VariantsList],
        max_neighbor_distance: int = MAX_NEIGHBOR_DISTANCE,
        num_threads: int = NUM_THREADS,
        match_all_breakpoints: bool = MATCH_ALL_BREAKPOINTS,
        match_variant_types: bool = MATCH_VARIANT_TYPES,
        min_ins_size_overlap: float = MIN_INS_SIZE_OVERLAP,
        min_del_size_overlap: float = MIN_DEL_SIZE_OVERLAP
) -> VariantsList:
    """
    Subtract variants lists from a target variants list.

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
        min_ins_size_overlap    :   Minimum insertion size overlap.
        min_del_size_overlap    :   Minimum deletion size overlap.

    Returns:
        VariantsList (with VariantCalls private to target_variants_list)
    """
    vl_subtracted = copy.deepcopy(target_variants_list)
    for variants_list in query_variants_lists:
        vl_subtracted = vl_subtracted.subtract(
            variants_list=variants_list,
            num_threads=num_threads,
            max_neighbor_distance=max_neighbor_distance,
            match_all_breakpoints=match_all_breakpoints,
            match_variant_types=match_variant_types,
            min_ins_size_overlap=min_ins_size_overlap,
            min_del_size_overlap=min_del_size_overlap
        )
    logger.info('%i variants and %i variant calls in the target VariantsList before diff.' %
                (len(target_variants_list.variants), len(target_variants_list.variant_call_ids)))
    logger.info('%i variants and %i variant calls in the target VariantsList after diff.' %
                (len(vl_subtracted.variants), len(vl_subtracted.variant_call_ids)))
    return vl_subtracted


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
    if variant_calling_method == VariantCallingMethods.CLAIRS:
        variants_list = parse_clairs_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id,
            case_id=case_id
        )
    elif variant_calling_method == VariantCallingMethods.CUTESV:
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
    elif variant_calling_method == VariantCallingMethods.MANTA_SOMATIC:
        variants_list = parse_manta_somatic_callset(
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
    elif variant_calling_method == VariantCallingMethods.SAVANA:
        variants_list = parse_savana_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id,
            case_id=case_id,
            control_id=control_id
        )
    elif variant_calling_method == VariantCallingMethods.SEVERUS:
        variants_list = parse_severus_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id,
            case_id=case_id
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
    elif variant_calling_method == VariantCallingMethods.SVISIONPRO:
        variants_list = parse_svisionpro_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id,
            case_id=case_id,
            control_id=control_id
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
