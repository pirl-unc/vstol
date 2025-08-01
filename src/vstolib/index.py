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
from intervaltree import IntervalTree
from typing import Dict


def index_genomic_ranges_list(df: pd.DataFrame, buffer: int) -> Dict[str, IntervalTree]:
    tree = {}
    for i in range(0, len(df)):
        row = df.iloc[i]
        if row['chromosome'] not in tree:
            tree[row['chromosome']] = IntervalTree()
        start = row['start'] - buffer
        end = row['end'] + buffer
        tree[row['chromosome']][start:end] = i

    return tree


def index_variants_list(df: pd.DataFrame, buffer: int) -> Dict[str, IntervalTree]:
    tree = {}
    for i in range(0, len(df)):
        row = df.iloc[i]
        if row['chromosome_1'] not in tree:
            tree[row['chromosome_1']] = IntervalTree()
        start = row['position_1'] - buffer
        end = row['position_1'] + buffer
        tree[row['chromosome_1']][start:end] = i

        if row['chromosome_2'] not in tree:
            tree[row['chromosome_2']] = IntervalTree()
        start = row['position_2'] - buffer
        end = row['position_2'] + buffer
        tree[row['chromosome_2']][start:end] = i

    return tree
