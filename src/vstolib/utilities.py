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


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise Exception('Boolean value expected.')

