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
The purpose of this python3 script is to implement the Annotator dataclass.
"""


import pandas as pd
from dataclasses import dataclass
from .logging import get_logger


logger = get_logger(__name__)


@dataclass
class Annotator:

    def annotate(self, df_variants: pd.DataFrame, num_processes: int) -> pd.DataFrame:
        raise Exception("Subclass must implement 'annotate' method")

