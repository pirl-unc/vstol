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
from typing import Optional
from networkx.utils.union_find import UnionFind
from .annotator import Annotator
from .collapse import collapse_group
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


def collapse(
        df_variants: pd.DataFrame,
        tumor_sample_ids: List[str],
        normal_sample_ids: Optional[List[str]] = None,
        min_methods: int = MIN_METHODS
) -> pd.DataFrame:
    """
    Collapse a merged variants Pandas DataFrame into a single representative
    row per `group_id`, applying somatic filters along the way.

    Filters (in order):
        1. (Only if `normal_sample_ids` is non-empty.) Drop any group in
           which at least one row whose `sample_id` is in
           `normal_sample_ids` has `position_1_read_count_alternate_allele > 0`
           or `position_2_read_count_alternate_allele > 0` (normal contamination).
        2. Drop any group whose number of unique `method` values is below
           `min_methods`.

    Collapse strategy:
        - Variant identity columns (chromosome_*, strand_*, operation_*,
          variant_type, sequence, filter) are taken from the first tumor row in
          the group. In particular `sequence` and `filter` are kept from that
          representative row only — they are never pipe-joined.
        - `position_1`, `position_2`, and `variant_size` are the median
          across tumor rows (rounded to int) — robust to caller-specific
          breakpoint jitter. Each position's median is taken only over rows
          sharing the representative breakend orientation (same chromosome),
          so a chromosome is never paired with a position from another
          chromosome (e.g. for translocations grouped from both mate orders).
        - Numeric evidence columns (quality, read counts, AF) are the max
          across tumor rows — the most supportive caller's value.
        - All remaining columns are pipe-joined unique values across the
          whole group (dedup with first-occurrence order), so no original
          information is lost.
        - Two informational columns are appended: `num_methods` and
          `num_rows_in_group`.
        - Output column renames: each surviving group becomes one
          representative variant, so the `group_id` column is renamed to
          `variant_id`, and the original per-row `variant_id` values
          (pipe-joined across the group) are exposed as `variant_ids`.

    Parameters:
        df_variants         :   Pandas DataFrame with a `group_id` column
                                (output of `vstol merge`).
        tumor_sample_ids    :   Sample IDs treated as tumor. Used to pick
                                the canonical variant representation.
        normal_sample_ids   :   Sample IDs treated as normal. Used for the
                                contamination filter. May be None or empty
                                if the merged TSV contains no normal sample
                                (e.g. tumor-only callsets), in which case
                                the contamination filter is skipped.
        min_methods         :   Minimum number of unique calling methods
                                required per group.

    Returns:
        Pandas DataFrame with one row per surviving group.
    """
    if 'group_id' not in df_variants.columns:
        raise ValueError(
            "Input DataFrame has no 'group_id' column. "
            "`collapse` expects the output of `vstol merge`."
        )

    if normal_sample_ids is None:
        normal_sample_ids = []

    overlap = set(tumor_sample_ids) & set(normal_sample_ids)
    if overlap:
        raise ValueError(
            "Sample IDs appear in both --tumor-sample-ids and --normal-sample-ids: %s"
            % sorted(overlap)
        )

    num_groups_in = df_variants['group_id'].nunique()

    # Step 1. Normal contamination filter (skipped if no normal samples given)
    if normal_sample_ids:
        normal_mask = df_variants['sample_id'].isin(normal_sample_ids)
        normal_has_alt = normal_mask & (
            (df_variants['position_1_read_count_alternate_allele'] > 0) |
            (df_variants['position_2_read_count_alternate_allele'] > 0)
        )
        contaminated_group_ids = set(df_variants.loc[normal_has_alt, 'group_id'].unique())
        df_clean = df_variants[~df_variants['group_id'].isin(contaminated_group_ids)]
        logger.info(
            'Dropped %i/%i group(s) due to normal contamination.' %
            (len(contaminated_group_ids), num_groups_in)
        )
    else:
        df_clean = df_variants
        logger.info(
            'No normal sample IDs given; skipping normal contamination filter.'
        )

    # Step 2. Min methods filter
    methods_per_group = df_clean.groupby('group_id')['method'].nunique()
    valid_group_ids = set(methods_per_group[methods_per_group >= min_methods].index)
    dropped_by_min_methods = set(methods_per_group.index) - valid_group_ids
    df_filtered = df_clean[df_clean['group_id'].isin(valid_group_ids)]
    logger.info(
        'Dropped %i group(s) with fewer than %i unique method(s).' %
        (len(dropped_by_min_methods), min_methods)
    )

    # Step 3. Collapse each surviving group
    columns = list(df_variants.columns)
    collapsed_rows = []
    dropped_no_tumor = 0
    for group_id, df_group in df_filtered.groupby('group_id', sort=True):
        df_tumor = df_group[df_group['sample_id'].isin(tumor_sample_ids)]
        if df_tumor.empty:
            dropped_no_tumor += 1
            logger.warning(
                'Group %s has no tumor rows; dropping.' % str(group_id)
            )
            continue
        collapsed_rows.append(
            collapse_group(group_id, df_group, df_tumor, columns)
        )
    if dropped_no_tumor:
        logger.info(
            'Dropped %i group(s) with no tumor sample rows.' % dropped_no_tumor
        )
    logger.info(
        'Collapsed %i group(s) from %i input group(s).' %
        (len(collapsed_rows), num_groups_in)
    )

    if not collapsed_rows:
        # Preserve the expected column order even if the result is empty.
        empty_cols = columns + ['num_methods', 'num_rows_in_group']
        return pd.DataFrame({c: [] for c in empty_cols}).rename(
            columns={
                'variant_id': 'variant_ids',
                'group_id': 'variant_id'
            }
        )

    df_collapsed = pd.DataFrame(collapsed_rows)
    # Preserve input column order, then append the new informational columns.
    ordered_cols = [c for c in columns if c in df_collapsed.columns] + \
                   ['num_methods', 'num_rows_in_group']
    # Expose the group identifier as `variant_id` and the joined per-row
    # `variant_id` values as `variant_ids`
    return df_collapsed[ordered_cols].rename(columns={
        'variant_id': 'variant_ids',
        'group_id': 'variant_id'
    })


def diff(
        df_variants: pd.DataFrame,
        df_variants_list: List[pd.DataFrame],
        match_both_positions: bool,
        max_breakpoint_distance: int,
        max_small_var_size: int,
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
        max_breakpoint_distance     :   Maximum breakpoint distance (applied to SVs).
        max_small_var_size          :   Variants with `0 < size <= max_small_var_size` must
                                        match positions exactly (no slop).
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
                max_breakpoint_distance,
                max_small_var_size
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
        max_small_var_size: int,
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
        max_breakpoint_distance     :   Maximum breakpoint distance (applied to SVs).
        max_small_var_size          :   Variants with `0 < size <= max_small_var_size` must
                                        match positions exactly (no slop).
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
                    max_breakpoint_distance,
                    max_small_var_size
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
    # Use the ordered union of columns across ALL input lists, not just list[0]. Inputs from
    # different callers/versions can have different column sets; keying off list[0] alone
    # raised KeyError (extra column) or produced ragged lists -> ValueError (missing column).
    columns = list(dict.fromkeys([c for df in df_variants_list for c in df.columns]))
    # `group_id` is assigned fresh below. If any input already carries a
    # `group_id` column (e.g. a previously-merged or collapsed TSV), drop it
    # here — otherwise it sits in both the explicit prefix and `columns` and is
    # appended twice per row, yielding unequal-length arrays in `pd.DataFrame`.
    columns = [c for c in columns if c != 'group_id']
    data = {k: [] for k in ['group_id'] + columns}
    group_id = 1
    for group in groups:
        for identifier in group:
            i, j = map(int, identifier.split(":"))
            row_dict = df_variants_list[i].iloc[j].to_dict()
            data['group_id'].append(group_id)
            for k in columns:
                data[k].append(row_dict.get(k))
        group_id += 1

    return pd.DataFrame(data)


def merge(
        df_variants_list: List[pd.DataFrame],
        match_both_positions: bool,
        max_breakpoint_distance: int,
        max_small_var_size: int,
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
        max_breakpoint_distance     :   Maximum breakpoint distance (applied to SVs).
        max_small_var_size          :   Variants with `0 < size <= max_small_var_size` must
                                        match positions exactly (no slop).
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
                    max_breakpoint_distance,
                    max_small_var_size
                )
            )

    with Pool(processes=num_processes) as pool:
        results = pool.starmap(intersect_variant_lists, args)

    # Step 2. Union-find from parallel outputs
    uf = UnionFind()
    for match_list in results:
        for id1, id2 in match_list:
            uf.union(id1, id2)

    # Step 2b. Union rows within the same dataframe that share the same `id`.
    # Some parsers (e.g. SAVANA) emit one row per sample for the same variant
    # call, all sharing a single `id` value, so those rows must end up in the
    # same group. This is O(M) per dataframe via groupby and avoids an
    # expensive self-pair `intersect_variant_lists` call.
    for i in range(0, len(df_variants_list)):
        for _, positions in df_variants_list[i].groupby('variant_id').indices.items():
            if len(positions) > 1:
                first = '%i:%i' % (i, int(positions[0]))
                for k in positions[1:]:
                    uf.union(first, '%i:%i' % (i, int(k)))

    # Step 3. Prepare the output DataFrame
    groups = uf.to_sets()
    # Ordered union of columns across ALL inputs (see note in intersect): different callers
    # can contribute different column sets, so keying off list[0] alone would crash.
    columns = list(dict.fromkeys([c for df in df_variants_list for c in df.columns]))
    # `group_id` is assigned fresh below. If any input already carries a
    # `group_id` column (e.g. a previously-merged or collapsed TSV), drop it
    # here — otherwise it sits in both the explicit prefix and `columns` and is
    # appended twice per row, yielding unequal-length arrays in `pd.DataFrame`.
    columns = [c for c in columns if c != 'group_id']
    data = {k: [] for k in ['group_id'] + columns}
    group_id = 1
    included_ids = set()
    for group in groups:
        for identifier in group:
            included_ids.add(identifier)
            i, j = map(int, identifier.split(":"))
            row_dict = df_variants_list[i].iloc[j].to_dict()
            data['group_id'].append(group_id)
            for k in columns:
                data[k].append(row_dict.get(k))
        group_id += 1
    for i in range(0, len(df_variants_list)):
        for j in range(0, len(df_variants_list[i])):
            if '%i:%i' % (i, j) not in included_ids:
                row_dict = df_variants_list[i].iloc[j].to_dict()
                data['group_id'].append(group_id)
                for k in columns:
                    data[k].append(row_dict.get(k))
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
    if method == VariantCallingMethod.CLAIR3:
        df_variants = parse_clair3_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.CLAIRS:
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
    elif method == VariantCallingMethod.DEEPSOMATIC:
        df_variants = parse_deepsomatic_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.DEEPVARIANT:
        df_variants = parse_deepvariant_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.DELLY2:
        df_variants = parse_delly2_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.DYSGU:
        df_variants = parse_dysgu_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.GRIDSS2:
        df_variants = parse_gridss2_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.HAPLOTYPECALLER:
        df_variants = parse_haplotypecaller_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.LONGSHOT:
        df_variants = parse_longshot_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.LUMPY:
        df_variants = parse_lumpy_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.MANTA:
        df_variants = parse_manta_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.MUTECT2:
        df_variants = parse_mutect2_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.NANOCALLER:
        df_variants = parse_nanocaller_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.NANOMONSV:
        df_variants = parse_nanomonsv_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.NANOVAR:
        df_variants = parse_nanovar_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.OCTOPUS:
        df_variants = parse_octopus_vcf(
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
    elif method == VariantCallingMethod.STRELKA2_GERMLINE:
        df_variants = parse_strelka2_germline_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.STRELKA2_SOMATIC:
        df_variants = parse_strelka2_somatic_vcf(
            df_vcf=df_vcf,
            platform=platform,
            source_id=source_id
        )
    elif method == VariantCallingMethod.SVABA:
        df_variants = parse_svaba_vcf(
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

