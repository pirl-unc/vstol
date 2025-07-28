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


from multiprocessing import Pool
from networkx.utils.union_find import UnionFind
from .annotator import Annotator
from .constants import VariantCallingMethod
from .default import *
from .intersect import *
from .intersect import intersect_variant_lists
from .logging import get_logger
from .overlap import overlap_variant_list
from .vcf import *


logger = get_logger(__name__)


def annotate(
        df_variants: pd.DataFrame,
        annotator: Annotator,
        num_processes: int
) -> pd.DataFrame:
    """
    Annotate a Pandas DataFrame of variants.

    Parameters:
        df_variants     :   Pandas DataFrame.
        annotator       :   Annotator.
        num_processes   :   Number of processes.

    Returns:
        Pandas DataFrame.
    """
    return annotator.annotate(
        df_variants=df_variants,
        num_processes=num_processes
    )


def diff(
        df_variants: pd.DataFrame,
        df_variants_list: List[pd.DataFrame],
        match_both_positions: bool,
        max_breakpoint_distance: int,
        match_operation_types: bool,
        num_processes: int = NUM_PROCESSES
) -> pd.DataFrame:
    """
    Diff a list of variants from a variant list.

    Parameters:
        df_variants                 :   Pandas DataFrame with the following columns:
                                        'chromosome_1',
                                        'position_1',
                                        'operation_1',
                                        'chromosome_2',
                                        'position_2',
                                        'operation_2'
        df_variants_list            :   List of Pandas DataFrames with the following columns:
                                        'chromosome_1',
                                        'position_1',
                                        'operation_1',
                                        'chromosome_2',
                                        'position_2',
                                        'operation_2'
        match_both_positions        :   If True, for two rows to be considered intersecting,
                                        both breakpoint positions must match or be near each other (max_breakpoint_distance).
        max_breakpoint_distance     :   Maximum breakpoint distance.
        match_operation_types       :   If True, operation_1 and operation_2 must match.
        num_processes               :   Number of processes.

    Returns:
        Pandas DataFrame with the following column appended:
        'group_id'
    """
    # Step 1. Identify intersections
    args = []
    for i in range(0, len(df_variants_list)):
        args.append(
            (
                0,
                i + 1,
                df_variants,
                df_variants_list[i],
                match_both_positions,
                match_operation_types,
                max_breakpoint_distance
            )
        )

    with Pool(processes=num_processes) as pool:
        results = pool.starmap(intersect_variant_lists, args)

    # Step 2. Identify dataframe indices to filter out
    positions_to_remove = set()
    for match_list in results:
        for id1, id2 in match_list:
            idx = int(id1.split(':')[1])
            positions_to_remove.add(idx)

    # Step 3. Prepare the output DataFrame
    df_variants_filtered = df_variants.drop(df_variants.index[list(positions_to_remove)])

    return df_variants_filtered


def intersect(
        df_variants_list: List[pd.DataFrame],
        match_both_positions: bool,
        max_breakpoint_distance: int,
        match_operation_types: bool,
        num_processes: int = NUM_PROCESSES
) -> pd.DataFrame:
    """
    Intersect variant lists.

    Parameters:
        df_variants_list            :   List of Pandas DataFrames with the following columns:
                                        'chromosome_1',
                                        'position_1',
                                        'operation_1',
                                        'chromosome_2',
                                        'position_2',
                                        'operation_2'
        match_both_positions        :   If True, for two rows to be considered intersecting,
                                        both breakpoint positions must match or be near each other (max_breakpoint_distance).
        max_breakpoint_distance     :   Maximum breakpoint distance.
        match_operation_types       :   If True, operation_1 and operation_2 must match.
        num_processes               :   Number of processes.

    Returns:
        Pandas DataFrame with the following column appended:
        'group_id'
    """
    # Step 1. Identify intersections
    args = []
    for i in range(len(df_variants_list)):
        for j in range(i + 1, len(df_variants_list)):
            args.append(
                (
                    i,
                    j,
                    df_variants_list[i],
                    df_variants_list[j],
                    match_both_positions,
                    match_operation_types,
                    max_breakpoint_distance
                )
            )

    with Pool(processes=num_processes) as pool:
        results = pool.starmap(intersect_variant_lists, args)

    # Step 2. Union-find from parallel outputs
    uf = UnionFind()
    for match_list in results:
        for id1, id2 in match_list:
            uf.union(id1, id2)

    # Step 3. Prepare the output DataFrame
    groups = uf.to_sets()
    data = {k: [] for k in ['group_id'] + list(df_variants_list[0].columns)}
    group_id = 1
    for group in groups:
        for identifier in group:
            i, j = map(int, identifier.split(":"))
            data['group_id'].append(group_id)
            for k, v in df_variants_list[i].iloc[j].to_dict().items():
                data[k].append(v)
        group_id += 1

    return pd.DataFrame(data)


def merge(
        df_variants_list: List[pd.DataFrame],
        match_both_positions: bool,
        max_breakpoint_distance: int,
        match_operation_types: bool,
        num_processes: int = NUM_PROCESSES
) -> pd.DataFrame:
    """
    Merge variant lists.

    Parameters:
        df_variants_list            :   List of Pandas DataFrames with the following columns:
                                        'chromosome_1',
                                        'position_1',
                                        'chromosome_2',
                                        'position_2'
        match_both_positions        :   If True, for two rows to be considered intersecting,
                                        both breakpoint positions must match or be near each other (max_breakpoint_distance).
        max_breakpoint_distance     :   Maximum breakpoint distance.
        match_operation_types       :   If True, operation_1 and operation_2 must match.
        num_processes               :   Number of processes.

    Returns:
        Pandas DataFrame with the following column appended:
        'group_id'
    """
    # Step 1. Identify intersections
    args = []
    for i in range(0, len(df_variants_list)):
        for j in range(i + 1, len(df_variants_list)):
            args.append(
                (
                    i,
                    j,
                    df_variants_list[i],
                    df_variants_list[j],
                    match_both_positions,
                    match_operation_types,
                    max_breakpoint_distance
                )
            )

    with Pool(processes=num_processes) as pool:
        results = pool.starmap(intersect_variant_lists, args)

    # Step 2. Union-find from parallel outputs
    uf = UnionFind()
    for match_list in results:
        for id1, id2 in match_list:
            uf.union(id1, id2)

    # Step 3. Prepare the output DataFrame
    groups = uf.to_sets()
    data = {k: [] for k in ['group_id'] + list(df_variants_list[0].columns)}
    group_id = 1
    included_ids = set()
    for group in groups:
        for identifier in group:
            included_ids.add(identifier)
            i, j = map(int, identifier.split(":"))
            data['group_id'].append(group_id)
            for k, v in df_variants_list[i].iloc[j].to_dict().items():
                data[k].append(v)
        group_id += 1
    for i in range(0, len(df_variants_list)):
        for j in range(0, len(df_variants_list[i])):
            if '%i:%i' % (i, j) not in included_ids:
                data['group_id'].append(group_id)
                for k, v in df_variants_list[i].iloc[j].to_dict().items():
                    data[k].append(v)
                group_id += 1

    return pd.DataFrame(data)


def overlap(
        df_variants: pd.DataFrame,
        df_genomic_ranges: pd.DataFrame,
        buffer: int
) -> pd.DataFrame:
    """
    Overlap variant list.

    Parameters:
        df_variant_list     :   Pandas DataFrame with the following columns:
                                'chromosome_1',
                                'position_1',
                                'chromosome_2',
                                'position_2'
        df_genomic_ranges   :   Pandas DataFrame with the following columns:
                                'chromosome',
                                'start',
                                'end'
        buffer              :   Buffer.

    Returns:
        Pandas DataFrame of variants overlapping any of the genomic ranges.
    """
    return overlap_variant_list(
        df_variants=df_variants,
        df_genomic_ranges=df_genomic_ranges,
        buffer=buffer
    )


def vcf2tsv(
        df_vcf: pd.DataFrame,
        source_id: str,
        method: str,
        platform: str
) -> pd.DataFrame:
    """
    Convert a VCF Pandas DataFrame to a VariantsList object.

    Parameters:
        df_vcf                  :   Pandas DataFrame (read from vcf.common.read_vcf_file).
        source_id               :   Source ID (e.g. patient ID or cell line sample ID).
        method                  :   Variant calling method.
        platform                :   Sequencing platform.

    Returns:
        Pandas DataFrame
    """
    if method == VariantCallingMethod.CLAIRS:
        df_variants = parse_clairs_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.CUTESV:
        df_variants = parse_cutesv_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.DBSNP:
        df_variants = parse_dbsnp_vcf(
            df_vcf=df_vcf
        )
    elif method == VariantCallingMethod.DEEPVARIANT:
        df_variants = parse_deepvariant_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.DELLY2_SOMATIC:
        df_variants = parse_delly2_somatic_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.GATK4_MUTECT2:
        df_variants = parse_gatk4_mutect2_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.LUMPY_SOMATIC:
        df_variants = parse_lumpy_somatic_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.MANTA_SOMATIC:
        df_variants = parse_manta_somatic_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.PBSV:
        df_variants = parse_pbsv_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.SAVANA:
        df_variants = parse_savana_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.SEVERUS:
        df_variants = parse_severus_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.SNIFFLES2:
        df_variants = parse_sniffles2_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.STRELKA2_SOMATIC:
        df_variants = parse_strelka2_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.SVIM:
        df_variants = parse_svim_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.SVISIONPRO:
        df_variants = parse_svisionpro_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    else:
        raise Exception('Unsupported variant calling method: %s' % method)
    return df_variants

