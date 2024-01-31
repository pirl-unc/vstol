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
The purpose of this python3 script is to define Velo default parameters.
"""


RANDOM_SEED = 1
NUM_THREADS = 4


"""annotate"""
# ANNOVAR protocol and corresponding operation
ANNOVAR_PROTOCOL_OPERATION = {
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


"""diff"""
DIFF_MAX_NEIGHBOR_DISTANCE = 1000
DIFF_MATCH_ALL_BREAKPOINTS = 'no'
DIFF_MATCH_VARIANT_TYPES = 'no'


"""filter"""
# Padding for an excluded region.
# The pad is applied to upstream and downstream of a gapped genomic region.
FILTER_EXCLUDED_REGION_PADDING = 100000

# Homopolymer length.
FILTER_HOMOPOLYMER_LENGTH = 20


"""intersect"""
# Maximum neighbor distance (bases).
INTERSECT_MAX_NEIGHBOR_DISTANCE = 10
INTERSECT_MATCH_ALL_BREAKPOINTS = 'yes'
INTERSECT_MATCH_VARIANT_TYPES = 'yes'



"""merge"""
# Maximum neighbor distance (bases).
MERGE_MAX_NEIGHBOR_DISTANCE = 10
MERGE_MATCH_ALL_BREAKPOINTS = 'yes'
MERGE_MATCH_VARIANT_TYPES = 'yes'


"""overlap"""
OVERLAP_PADDING = 0


"""tsv2vcf"""


"""vcf2tsv"""



