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
The purpose of this python3 script is to define VSTOL default parameters.
"""


ANNOVAR_PROTOCOL_OPERATION = {     # ANNOVAR protocol and corresponding operation
    'refGene': 'g',
    'exac03': 'f',
    '1000g2015aug_eur': 'f',
    '1000g2015aug_eas': 'f',
    '1000g2015aug_sas': 'f',
    'clinvar_20210501': 'f',
    'cosmic96_coding': 'f',
    'avsnp150': 'f',
    'dbnsfp42c': 'f'
}
GZIP = 'no'
HOMOPOLYMER_LENGTH = 10
MATCH_ALL_BREAKPOINTS = 'yes'
MATCH_VARIANT_TYPES = 'yes'
MAX_NEIGHBOR_DISTANCE = 100
MIN_DEL_SIZE_OVERLAP = 0.5
MIN_INS_SIZE_OVERLAP = 0.5
NUM_THREADS = 4
RANDOM_SEED = 1
WINDOW = 10000



