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
from typing import  List, Tuple
from .index import index_variants_list


def intersect_variant_lists(
    i: int,
    j: int,
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    match_both_positions: bool,
    match_operation_types: bool,
    max_breakpoint_distance: int
) -> List[Tuple[str, str]]:
    """
    Find intersections.

    Parameters:
        i                           :   DataFrame 1 ID.
        j                           :   DataFrame 2 ID.
        df1                         :   Pandas DataFrame 1.
        df2                         :   Pandas DataFrame 2.
        match_both_positions        :   If True, for two rows to be considered intersecting,
                                        both breakpoint positions must match or be near each other (max_breakpoint_distance).
        max_breakpoint_distance     :   Maximum breakpoint distance.
        match_operation_types       :   If True, operation_1 and operation_2 must match.

    Returns:
        List of intersecting pairs (<dataframe_id>:<dataframe_position>, <dataframe_id>:<dataframe_position>)
    """
    tree_dict = index_variants_list(df=df2, buffer=max_breakpoint_distance)
    matches = []
    for index in range(len(df1)):
        row = df1.iloc[index]
        if row['chromosome_1'] in tree_dict.keys():
            matches_1 = tree_dict[row['chromosome_1']][row['position_1']]
        else:
            matches_1 = set([])
        if row['chromosome_2'] in tree_dict.keys():
            matches_2 = tree_dict[row['chromosome_2']][row['position_2']]
        else:
            matches_2 = set([])
        for match in matches_1.union(matches_2):
            row2 = df2.iloc[match.data]
            id1 = f"{i}:{index}"
            id2 = f"{j}:{match.data}"
            if match_both_positions and match_operation_types:
                intersects = (
                     (row2['chromosome_1'] == row['chromosome_1']) &
                     (row2['position_1'] >= row['position_1'] - max_breakpoint_distance) &
                     (row2['position_1'] <= row['position_1'] + max_breakpoint_distance) &
                     (row2['operation_1'] == row['operation_1']) &
                     (row2['chromosome_2'] == row['chromosome_2']) &
                     (row2['position_2'] <= row['position_2'] + max_breakpoint_distance) &
                     (row2['position_2'] >= row['position_2'] - max_breakpoint_distance) &
                     (row2['operation_2'] == row['operation_2'])
                ) or (
                     (row2['chromosome_1'] == row['chromosome_2']) &
                     (row2['position_1'] >= row['position_2'] - max_breakpoint_distance) &
                     (row2['position_1'] <= row['position_2'] + max_breakpoint_distance) &
                     (row2['operation_1'] == row['operation_2']) &
                     (row2['chromosome_2'] == row['chromosome_1']) &
                     (row2['position_2'] >= row['position_1'] - max_breakpoint_distance) &
                     (row2['position_2'] <= row['position_1'] + max_breakpoint_distance) &
                     (row2['operation_2'] == row['operation_1'])
                )
            elif match_both_positions == False and match_operation_types:
                intersects = (
                     (row2['chromosome_1'] == row['chromosome_1']) &
                     (row2['position_1'] >= row['position_1'] - max_breakpoint_distance) &
                     (row2['position_1'] <= row['position_1'] + max_breakpoint_distance) &
                     (row2['operation_1'] == row['operation_1'])
                ) or (
                     (row2['chromosome_2'] == row['chromosome_2']) &
                     (row2['position_2'] <= row['position_2'] + max_breakpoint_distance) &
                     (row2['position_2'] >= row['position_2'] - max_breakpoint_distance) &
                     (row2['operation_2'] == row['operation_2'])
                ) or (
                     (row2['chromosome_1'] == row['chromosome_2']) &
                     (row2['position_1'] >= row['position_2'] - max_breakpoint_distance) &
                     (row2['position_1'] <= row['position_2'] + max_breakpoint_distance) &
                     (row2['operation_1'] == row['operation_2'])
                ) or (
                     (row2['chromosome_2'] == row['chromosome_1']) &
                     (row2['position_2'] >= row['position_1'] - max_breakpoint_distance) &
                     (row2['position_2'] <= row['position_1'] + max_breakpoint_distance) &
                     (row2['operation_2'] == row['operation_1'])
                )
            elif match_both_positions and match_operation_types == False:
                intersects = (
                     (row2['chromosome_1'] == row['chromosome_1']) &
                     (row2['position_1'] >= row['position_1'] - max_breakpoint_distance) &
                     (row2['position_1'] <= row['position_1'] + max_breakpoint_distance) &
                     (row2['chromosome_2'] == row['chromosome_2']) &
                     (row2['position_2'] <= row['position_2'] + max_breakpoint_distance) &
                     (row2['position_2'] >= row['position_2'] - max_breakpoint_distance)
                ) or (
                     (row2['chromosome_1'] == row['chromosome_2']) &
                     (row2['position_1'] >= row['position_2'] - max_breakpoint_distance) &
                     (row2['position_1'] <= row['position_2'] + max_breakpoint_distance) &
                     (row2['chromosome_2'] == row['chromosome_1']) &
                     (row2['position_2'] >= row['position_1'] - max_breakpoint_distance) &
                     (row2['position_2'] <= row['position_1'] + max_breakpoint_distance)
                )
            else:
                intersects = (
                     (row2['chromosome_1'] == row['chromosome_1']) &
                     (row2['position_1'] >= row['position_1'] - max_breakpoint_distance) &
                     (row2['position_1'] <= row['position_1'] + max_breakpoint_distance)
                ) or (
                     (row2['chromosome_2'] == row['chromosome_2']) &
                     (row2['position_2'] <= row['position_2'] + max_breakpoint_distance) &
                     (row2['position_2'] >= row['position_2'] - max_breakpoint_distance)
                ) or (
                     (row2['chromosome_1'] == row['chromosome_2']) &
                     (row2['position_1'] >= row['position_2'] - max_breakpoint_distance) &
                     (row2['position_1'] <= row['position_2'] + max_breakpoint_distance)
                ) or (
                     (row2['chromosome_2'] == row['chromosome_1']) &
                     (row2['position_2'] >= row['position_1'] - max_breakpoint_distance) &
                     (row2['position_2'] <= row['position_1'] + max_breakpoint_distance)
                )
            if intersects:
                matches.append((id1, id2))

    return matches

