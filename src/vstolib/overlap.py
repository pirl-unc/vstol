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
from .index import index_genomic_ranges_list


def overlap_variant_list(
    df_variants: pd.DataFrame,
    df_genomic_ranges: pd.DataFrame,
    buffer: int
) -> pd.DataFrame:
    """
    Overlap variant list.

    Parameters:
        df_variants         :   Pandas DataFrame of variants.
        df_genomic_ranges   :   Pandas DataFrame of genomic ranges.
        buffer              :   Buffer.

    Returns:
        Pandas DataFrame of overlapped variants.
    """
    tree_dict = index_genomic_ranges_list(df=df_genomic_ranges, buffer=buffer)
    overlapping_indices = set()
    for i in range(0, len(df_variants)):
        row = df_variants.iloc[i]
        # Guard chromosome membership: a variant on a contig with no genomic range
        # (the ranges file need not cover every contig) previously raised KeyError and
        # aborted the whole command. Treat a missing chromosome as zero overlaps.
        if row['chromosome_1'] in tree_dict:
            matches_1 = tree_dict[row['chromosome_1']][row['position_1']]
        else:
            matches_1 = set()
        if row['chromosome_2'] in tree_dict:
            matches_2 = tree_dict[row['chromosome_2']][row['position_2']]
        else:
            matches_2 = set()
        if len(matches_1) > 0 or len(matches_2) > 0:
            overlapping_indices.add(i)

    return df_variants.iloc[list(overlapping_indices)]

