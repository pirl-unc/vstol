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


from dataclasses import dataclass
from .logging import get_logger
from .variants_list import VariantsList


logger = get_logger(__name__)


@dataclass
class Annotator:

    def annotate(self, variants_list: VariantsList, num_processes: int) -> VariantsList:
        raise Exception("Subclass must implement 'annotate' method")

