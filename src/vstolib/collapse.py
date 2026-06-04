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


import pandas as pd
from typing import List


# Per-column collapse rules for `collapse()`. Any column not listed here
# defaults to 'join' (pipe-joined unique values, preserving first-occurrence
# order) so new parser-emitted columns degrade gracefully.
_COLLAPSE_RULES = {
    # Variant identity — first value from a representative (tumor) row.
    'chromosome_1': 'first',
    'chromosome_2': 'first',
    'strand_1': 'first',
    'strand_2': 'first',
    'operation_1': 'first',
    'operation_2': 'first',
    'variant_type': 'first',
    # Inserted/breakend sequence — keep the representative row's value only;
    # never pipe-join (a multi-sequence cell breaks downstream loaders, and a
    # concatenation of sequences from different rows is not itself a valid
    # sequence).
    'sequence': 'first',
    # Filter / QC status — keep the representative row's value only; never
    # pipe-join (a cell like "PASS|." is not a meaningful filter value and
    # breaks downstream loaders that expect a single token).
    'filter': 'first',
    # Positional / size — median across tumor rows, rounded to int.
    'position_1': 'median_int',
    'position_2': 'median_int',
    'variant_size': 'median_int',
    # Numeric evidence — max across tumor rows.
    'quality': 'max',
    'position_1_read_count_reference_allele': 'max',
    'position_2_read_count_reference_allele': 'max',
    'position_1_read_count_alternate_allele': 'max',
    'position_2_read_count_alternate_allele': 'max',
    'position_1_read_count_total': 'max',
    'position_2_read_count_total': 'max',
    'position_1_alternate_allele_fraction': 'max',
    'position_2_alternate_allele_fraction': 'max',
}


def _join_unique(series: pd.Series) -> str:
    """Pipe-join unique non-null values in first-occurrence order.

    Nulls are dropped first: a column that is entirely missing for a group
    (e.g. `sequence` or `read_ids_alternate_allele` for structural variants)
    collapses to an empty string rather than raising, and partially-missing
    columns never emit a literal ``nan`` token.
    """
    return '|'.join(series.dropna().astype(str).drop_duplicates().tolist())


# Each median-aggregated breakpoint position is only meaningful alongside rows
# sharing its chromosome. A group can hold the same event in opposite breakend
# orders (e.g. translocation mates chrA<->chrB and chrB<->chrA, or different
# callers ordering breakends differently); medianing `position_1` across all of
# them would average a chrA coordinate with a chrB one, leaving `chromosome_1`
# (taken from the representative row) paired with a position on neither
# chromosome. So restrict each position's median to rows whose governing
# chromosome matches the representative row's.
_POSITION_CHROMOSOME = {
    'position_1': 'chromosome_1',
    'position_2': 'chromosome_2'
}


def collapse_group(
        group_id,
        df_group: pd.DataFrame,
        df_tumor: pd.DataFrame,
        columns: List[str]
) -> dict:
    """Collapse a single group into one representative row."""
    row = {'group_id': group_id}
    # The first tumor row defines the representative breakend orientation; the
    # 'first' identity columns and the position medians below are all keyed off
    # it so position stays paired with its own chromosome.
    rep = df_tumor.iloc[0]
    for col in columns:
        if col == 'group_id':
            continue
        rule = _COLLAPSE_RULES.get(col, 'join')
        if rule == 'first':
            row[col] = rep[col]
        elif rule == 'median_int':
            chrom_col = _POSITION_CHROMOSOME.get(col)
            if chrom_col is not None:
                # Median only over rows on the representative's chromosome.
                series = df_tumor.loc[df_tumor[chrom_col] == rep[chrom_col], col]
            else:
                series = df_tumor[col]
            row[col] = int(round(series.median()))
        elif rule == 'max':
            row[col] = df_tumor[col].max()
        else:  # 'join' — preserve all original values across the whole group.
            row[col] = _join_unique(df_group[col])
    row['num_methods'] = df_group['method'].nunique()
    row['num_rows_in_group'] = len(df_group)
    return row