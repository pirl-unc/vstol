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
The purpose of this python3 script is to implement general-purpose utility functions.
"""


import pandas as pd
from typing import Any, Dict, Literal
from .logging import get_logger


logger = get_logger(__name__)


def get_typed_value(
        value: Any,
        default_value: Any,
        type: Literal[int,float,str,bool]
) -> Any:
    """
    Safely converts a value from a VCF row.

    Parameters:
        value           :   Value.
        default_value   :   Default value.
        type            :   Type (str, int, float, or bool).

    Returns:
        Any
    """
    try:
        if pd.isna(value):
            value = default_value
        else:
            if type == str:
                value = str(value)
            elif type == int:
                value = int(value)
            elif type == float:
                value = float(value)
            elif type == bool:
                value = bool(value)
            else:
                value = default_value
    except:
        value = default_value
    return value


def is_integer(x) -> bool:
    try:
        int(x)
        return True
    except ValueError:
        return False
    except TypeError:
        return False


def is_repeated_sequence(sequence: str) -> bool:
    if sequence == '':
        return False
    sequence = sequence.upper()
    first_char = sequence[0]
    for char in sequence[1:]:
        if char != first_char:
            return False
    return True


def overlaps(
        df: pd.DataFrame,
        chrom: str,
        start: int,
        end: int
    ) -> bool:
    """
    Returns true if a genomic region overlaps any regions in a DataFrame.

    Parameters:
        df          :   DataFrame with the following columns:
                        'chromosome_1', 'position_1', 'chromosome_2', 'position_2'
        chrom       :   Chromosome.
        start       :   Start position.
        end         :   End position.

    Returns:
        True if any row in the DataFrame overlaps the queried region.
        False otherwise.
    """
    # De Morgan's law on checking for non-overlapping regions
    df_matched = df.loc[
        (df['chr_1'] == chrom) &
        (df['chr_2'] == chrom) &
        (df['pos_2'] >= start) &
        (df['pos_1'] <= end),:
    ]
    if len(df_matched) > 0:
        return True
    else:
        return False


def retrieve_from_dict(
        dct: Dict,
        key: str,
        default_value: Any,
        type: Literal[int,float,str,bool]
) -> Any:
    """
    Safely retrieves a value from a dictionary.

    Parameters:
        dct             :   Dictionary (or vector).
        key             :   Key.
        default_value   :   Default value.
        type            :   Type (int, float, str, bool).

    Returns:
        Any
    """
    try:
        value = dct[key]
    except:
        value = default_value
    return get_typed_value(value=value,
                           default_value=default_value,
                           type=type)


def retrieve_from_list(
        lst: list,
        index: int,
        default_value: Any,
        type: Literal[int,float,str,bool]
) -> Any:
    """
    Safely retrieves a value from a list.

    Parameters:
        lst             :   List.
        index           :   Key.
        default_value   :   Default value.
        type            :   Type (int, float, str, bool).

    Returns:
        Any
    """
    try:
        value = lst[index]
    except:
        value = default_value
    return get_typed_value(value=value,
                           default_value=default_value,
                           type=type)


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise Exception('Boolean value expected.')

