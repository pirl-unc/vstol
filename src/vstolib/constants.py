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
The purpose of this python3 script is to define constants.
"""


class Annotator:
    GENCODE = 'gencode'
    ALL = [
        GENCODE
    ]


class GenomicRegionTypes:
    EXONIC = 'exonic'
    INTRONIC = 'intronic'
    FIVE_PRIME_UTR = '5prime_utr'
    THREE_PRIME_UTR = '3prime_utr'
    UNTRANSLATED_REGION = 'utr'
    START_CODON = 'start_codon'
    STOP_CODON = 'stop_codon'
    INTERGENIC = 'intergenic'
    ALL = [
        EXONIC,
        INTRONIC,
        FIVE_PRIME_UTR,
        THREE_PRIME_UTR,
        UNTRANSLATED_REGION,
        START_CODON,
        STOP_CODON,
        INTERGENIC
    ]


class OperationType:
    UPSTREAM = 'U'
    DOWNSTREAM = 'D'
    UNKNOWN = ''
    ALL = [UPSTREAM,DOWNSTREAM,UNKNOWN]


class Strand:
    FORWARD = '+'
    REVERSE = '-'
    BOTH = '*'
    UNKNOWN = ''
    ALL = [FORWARD,REVERSE,BOTH,UNKNOWN]


class VariantCallingMethod:
    CLAIR3 = 'clair3'
    CLAIRS = 'clairs'
    CUTESV = 'cutesv'
    DBSNP = 'dbsnp'
    DEEPSOMATIC = 'deepsomatic'
    DEEPVARIANT = 'deepvariant'
    DELLY2 = 'delly2'
    DYSGU = 'dysgu'
    GRIDSS2 = 'gridss2'
    HAPLOTYPECALLER = 'haplotypecaller'
    LONGSHOT = 'longshot'
    LUMPY = 'lumpy'
    MANTA = 'manta'
    MUTECT2 = 'mutect2'
    NANOCALLER = 'nanocaller'
    NANOMONSV = 'nanomonsv'
    NANOVAR = 'nanovar'
    PBSV = 'pbsv'
    OCTOPUS = 'octopus'
    SAVANA = 'savana'
    SEVERUS = 'severus'
    SNIFFLES2 = 'sniffles2'
    STRELKA2_GERMLINE = 'strelka2-germline'
    STRELKA2_SOMATIC = 'strelka2-somatic'
    SVABA = 'svaba'
    SVIM = 'svim'
    SVISIONPRO = 'svisionpro'
    ALL = [
        CLAIR3,
        CLAIRS,
        CUTESV,
        DBSNP,
        DEEPSOMATIC,
        DEEPVARIANT,
        DELLY2,
        DYSGU,
        GRIDSS2,
        HAPLOTYPECALLER,
        LONGSHOT,
        LUMPY,
        MANTA,
        MUTECT2,
        NANOCALLER,
        NANOMONSV,
        NANOVAR,
        OCTOPUS,
        PBSV,
        SAVANA,
        SEVERUS,
        SNIFFLES2,
        STRELKA2_GERMLINE,
        STRELKA2_SOMATIC,
        SVABA,
        SVIM,
        SVISIONPRO
    ]


class VariantType:
    SINGLE_NUCLEOTIDE_VARIANT = 'SNV'
    MULTI_NUCLEOTIDE_VARIANT = 'MNV'
    INSERTION = 'INS'
    DELETION = 'DEL'
    DUPLICATION = 'DUP'
    INVERSION = 'INV'
    TRANSLOCATION = 'TRA'
    BREAKPOINT = 'BND'
    REFERENCE = 'REF'

    ALL = [
        SINGLE_NUCLEOTIDE_VARIANT,
        MULTI_NUCLEOTIDE_VARIANT,
        INSERTION,
        DELETION,
        DUPLICATION,
        INVERSION,
        TRANSLOCATION,
        BREAKPOINT
    ]
