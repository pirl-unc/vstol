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


class Annotators:
    ENSEMBL = 'ensembl'
    GENCODE = 'gencode'
    REFSEQ = 'refseq'
    ANNOVAR = 'annovar'
    ALL = [
        ENSEMBL,
        GENCODE,
        REFSEQ,
        ANNOVAR
    ]


class CollapseStrategies:
    MAX_ALTERNATE_ALLELE_READ_COUNT = 'max_alternate_allele_read_count'
    ALL = [
        MAX_ALTERNATE_ALLELE_READ_COUNT
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


class NucleicAcidTypes:
    DNA = 'dna'
    RNA = 'rna'
    DNA_RNA = 'dna_rna'
    ALL = [
        DNA,
        RNA,
        DNA_RNA
    ]


class Strands:
    POSITIVE = '+'
    NEGATIVE = '-'
    BOTH_STRANDS = '+-'


class TranslocationOrientations:
    ORIENTATION_1 = 't[p['  # piece extending to the right of p is joined after t
    ORIENTATION_2 = 't]p]'  # reverse complement piece extending left of p is joined after t
    ORIENTATION_3 = ']p]t'  # piece extending to the left of p is joined before t
    ORIENTATION_4 = '[p[t'  # reverse complement extending right of p is joined before t
    ALL = [
        ORIENTATION_1,
        ORIENTATION_2,
        ORIENTATION_3,
        ORIENTATION_4
    ]


class VariantCallingMethods:
    CLAIRS = 'clairs'
    CUTESV = 'cutesv'
    DBSNP = 'dbsnp'
    DEEPVARIANT = 'deepvariant'
    DELLY2_SOMATIC = 'delly2-somatic'
    GATK4_MUTECT2 = 'gatk4-mutect2'
    LUMPY_SOMATIC = 'lumpy-somatic'
    MANTA_SOMATIC = 'manta-somatic'
    PBSV = 'pbsv'
    SAVANA = 'savana'
    SEVERUS = 'severus'
    SNIFFLES2 = 'sniffles2'
    STRELKA2_SOMATIC = 'strelka2-somatic'
    SVIM = 'svim'
    SVISIONPRO = 'svisionpro'
    ALL = [
        CLAIRS,
        CUTESV,
        DBSNP,
        DEEPVARIANT,
        DELLY2_SOMATIC,
        GATK4_MUTECT2,
        LUMPY_SOMATIC,
        MANTA_SOMATIC,
        PBSV,
        SAVANA,
        SEVERUS,
        SNIFFLES2,
        STRELKA2_SOMATIC,
        SVIM,
        SVISIONPRO
    ]

    class AttributeTypes:
        CLAIRS = {
            'ID': str,
            'H': bool,
            'FAU': int,
            'FCU': int,
            'FGU': int,
            'FTU': int,
            'RAU': int,
            'RCU': int,
            'RGU': int,
            'RTU': int,
            'GT': str,
            'GQ': int,
            'DP': int,
            'AF': float,
            'AD': str,
            'NAF': float,
            'NDP': int,
            'NAD': int,
            'AU': int,
            'CU': int,
            'GU': int,
            'TU': int,
            'NAU': int,
            'NCU': int,
            'NGU': int,
            'NTU': int
        }
        CUTESV = {
            'ID': str,
            'SVTYPE': str,
            'SVLEN': int,
            'CHR2': str,
            'END': int,
            'CIPOS': str,
            'CILEN': str,
            'RE': int,
            'STRAND': str,
            'RNAMES': str,
            'AF': float,
            'PRECISE': bool,
            'GT': str,
            'GQ': float,
            'PL': str,
            'DR': int,
            'DV': int
        }
        DBSNP = {
            'ID': str
        }
        DELLY2_SOMATIC = {
            'ID': str,
            'SVTYPE': str,
            'SVMETHOD': str,
            'SVLEN': int,
            'END': int,
            'CHR2': str,
            'POS2': int,
            'PE': int,
            'MAPQ': int,
            'CT': str,
            'CIPOS': str,
            'CIEND': str,
            'SRMAPQ': int,
            'INSLEN': int,
            'HOMLEN': int,
            'SR': int,
            'SRQ': int,
            'CONSENSUS': str,
            'CE': float,
            'CONSBP': int,
            'RDRATIO': float,
            'GT': str,
            'GL': str,
            'GQ': int,
            'FT': str,
            'RC': int,
            'RCL': int,
            'RCR': int,
            'RDCN': int,
            'DR': int,
            'DV': int,
            'RR': int,
            'RV': int,
            'PRECISE': bool,
            'SOMATIC': bool
        }
        DEEPVARIANT = {
            'ID': str,
            'END': int,
            'GT': str,
            'GQ': int,
            'DP': int,
            'MIN_DP': int,
            'AD': str,
            'VAF': float,
            'PL': str,
            'MED_DP': int
        }
        GATK4_MUTECT2 = {
            'ID': str,
            'AS_FilterStatus': str,
            'AS_SB_TABLE': str,
            'AS_UNIQ_ALT_READ_COUNT': int,
            'CONTQ': float,
            'ECNT': int,
            'GERMQ': int,
            'MBQ': int,
            'MFRL': int,
            'MMQ': int,
            'MPOS': int,
            'NALOD': float,
            'NCOUNT': int,
            'NLOD': float,
            'OCM': int,
            'PON': bool,
            'POPAF': float,
            'AF': float,
            'ROQ': float,
            'RPA': int,
            'RU': str,
            'SEQQ': int,
            'STR': bool,
            'STRANDQ': int,
            'STRQ': int,
            'TLOD': float,
            'AD': str,
            'DP': int,
            'F1R2': str,
            'F2R1': str,
            'FAD': str,
            'GQ': float,
            'GT': str,
            'PGT': str,
            'PID': str,
            'PL': int,
            'PS': int,
            'SB': str
        }
        GRIDSS = {
            'GT': str,
            'AF': float,
            'ANRP': int,
            'ANRPQ': float,
            'ANSR': int,
            'ANSRQ': float,
            'ASQ': float,
            'ASRP': int,
            'ASSR': int,
            'BAQ': float,
            'BASRP': int,
            'BASSR': int,
            'BQ': float,
            'BSC': int,
            'BSCQ': float,
            'BUM': int,
            'BUMQ': float,
            'BVF': int,
            'CASQ': float,
            'IC': int,
            'IQ': float,
            'QUAL': float,
            'RASQ': float,
            'REF': int,
            'REFPAIR': int,
            'RP': int,
            'RPQ': float,
            'SR': int,
            'SRQ': float,
            'VF': int,
            'ANRP': int,
            'ANRPQ': float,
            'ANSR': int,
            'ANSRQ': float,
            'AS': int,
            'ASC': str,
            'ASQ': float,
            'ASRP': int,
            'ASSR': int,
            'BA': int,
            'BAQ': float,
            'BASRP': int,
            'BASSR': int,
            'BEALN': str,
            'BEID': str,
            'BEIDH': int,
            'BEIDL': int,
            'BENAMES': str,
            'BMQ': float,
            'BMQN': float,
            'BMQX': float,
            'BPNAMES': str,
            'BQ': float,
            'BSC': int,
            'BSCQ': float,
            'BUM': int,
            'BUMQ': float,
            'BVF': int,
            'CAS': int,
            'CASQ': float,
            'CIPOS': str,
            'CIRPOS': str,
            'CQ': float,
            'EVENT': str,
            'HOMLEN': int,
            'HOMSEQ': str,
            'IC': int,
            'IHOMPOS': str,
            'INSRM': str,
            'INSRMP': float,
            'INSRMRC': str,
            'INSRMRO': str,
            'INSRMRT': str,
            'INSTAXID': int,
            'IQ': float,
            'MATEID': str,
            'MQ': float,
            'MQN': float,
            'MQX': float,
            'RAS': int,
            'RASQ': float,
            'REF': int,
            'REFPAIR': int,
            'RF': int,
            'RP': int,
            'RPQ': float,
            'RSI': int,
            'SB': float,
            'SC': str,
            'SI': str,
            'SR': int,
            'SRQ': float,
            'SVTYPE': str,
            'VF': int,
            'TAF': float,
            'LOCAL_LINKED_BY': str,
            'REMOTE_LINKED_BY': str
        }
        LUMPY_SOMATIC = {
            'ID': str,
            'SVTYPE': str,
            'STRANDS': str,
            'SVLEN': int,
            'END': int,
            'CIPOS': str,
            'CIEND': str,
            'CIPOS95': str,
            'CIEND95': str,
            'SU': int,
            'PE': int,
            'SR': int,
            'GT': str,
            'BD': int,
            'MATEID': str,
            'EVENT': int,
            'EV': str,
            'PRPOS': str,
            'PREND': str,
            'PRECISE': bool,
            'IMPRECISE': bool,
            'SECONDARY': bool
        }
        MANTA_SOMATIC = {
            'ID': str,
            'SVTYPE': str,
            'SVLEN': int,
            'END': int,
            'CIPOS': str,
            'CIEND': str,
            'CIGAR': str,
            'MATEID': str,
            'EVENT': str,
            'HOMLEN': int,
            'HOMSEQ': str,
            'SVINSLEN': int,
            'SVINSSEQ': str,
            'LEFT_SVINSSEQ': str,
            'RIGHT_SVINSSEQ': str,
            'BND_DEPTH': int,
            'MATE_BND_DEPTH': int,
            'PRECISE': bool,
            'IMPRECISE': bool,
            'SOMATIC': bool,
            'SOMATICSCORE': int,
            'JUNCTION_SOMATICSCORE': int,
            'PR': str,
            'SR': str
        }
        PBSV = {
            'ID': str,
            'SVTYPE': str,
            'END': int,
            'SVLEN': int,
            'SVANN': str,
            'CIPOS': str,
            'MATEID': str,
            'MATEDIST': int,
            'PRECISE': bool,
            'GT': str,
            'DP': int,
            'AD': str,
            'SAC': str,
            'NotFullySpanned': bool
        }
        SAVANA = {
            'ID': str,
            'SVTYPE': str,
            'MATEID': str,
            'NORMAL_SUPPORT': int,
            'TUMOUR_SUPPORT': int,
            'SVLEN': float,
            'BP_NOTATION': str,
            'ORIGINATING_CLUSTER': str,
            'END_CLUSTER': str,
            'ORIGIN_STARTS_STD_DEV': float,
            'ORIGIN_MAPQ_MEAN': float,
            'ORIGIN_EVENT_SIZE_STD_DEV': float,
            'ORIGIN_EVENT_SIZE_MEDIAN': float,
            'ORIGIN_EVENT_SIZE_MEAN': float,
            'END_STARTS_STD_DEV': float,
            'END_MAPQ_MEAN': float,
            'END_EVENT_SIZE_STD_DEV': float,
            'END_EVENT_SIZE_MEDIAN': float,
            'END_EVENT_SIZE_MEAN': float,
            'TUMOUR_DP': str,
            'NORMAL_DP': str,
            'GT': str,
            'PRECISE': bool,
            'CLASS': str
        }
        SEVERUS = {
            'ID': str,
            'SVTYPE': str,
            'SVLEN': int,
            'CHR2': str,
            'END': int,
            'STRANDS': str,
            'DETAILED_TYPE': str,
            'INSLEN': int,
            'MAPQ': int,
            'PHASESETID': str,
            'HP': int,
            'CLUSTERID': str,
            'INSSEQ': str,
            'MATE_ID': str,
            'INSIDE_VNTR': str,
            'ALINGED_POS': str,
            'GT': str,
            'GQ': int,
            'DR': int,
            'DV': int,
            'VAF': float,
            'hVAF': str,
            'PRECISE': bool
        }
        SNIFFLES2 = {
            'ID': str,
            'SVLEN': int,
            'SVTYPE': str,
            'CHR2': str,
            'SUPPORT': int,
            'SUPPORT_INLINE': int,
            'SUPPORT_LONG': int,
            'END': int,
            'STDEV_POS': float,
            'STDEV_LEN': float,
            'COVERAGE': str,
            'STRAND': str,
            'AC': int,
            'SUPP_VEC': str,
            'CONSENSUS_SUPPORT': int,
            'RNAMES': str,
            'AF': float,
            'NM': float,
            'PHASE': str,
            'GT': str,
            'GQ': int,
            'DR': int,
            'DV': int,
            'PRECISE': bool
        }
        STRELKA2_SOMATIC = {
            'ID': str,
            'QSS': int,
            'TQSS': int,
            'NT': str,
            'QSS_NT': int,
            'TQSS_NT': int,
            'SGT': str,
            'MQ': float,
            'MQ0': int,
            'ReadPosRankSum': float,
            'PNOISE': float,
            'PNOISE2': float,
            'SomaticEVS': float,
            'QSI': int,
            'TQSI': int,
            'QSI_NT': int,
            'TQSI_NT': int,
            'RU': str,
            'RC': int,
            'IC': int,
            'IHP': int,
            'SOMATIC': bool,
            'OVERLAP': bool,
            'FDP': int,
            'SDP': int,
            'SUBDP': int,
            'AU': str,
            'CU': str,
            'GU': str,
            'TU': str,
            'DP': int,
            'DP2': int,
            'TAR': str,
            'TIR': str,
            'TOR': str,
            'SNVSB': float,
            'DP50': float,
            'FDP50': float,
            'SUBDP50': float,
            'BCN50': float,
            'END': int,
            'SNVHPOL': int,
            'CIGAR': str,
            'REFREP': int,
            'IDREP': int,
            'BLOCKAVG_MIN30P3A': bool,
            'GT': str,
            'GQ': int,
            'GQX': int,
            'DPF': int,
            'MIN_DP': int,
            'AD': str,
            'ADF': str,
            'ADR': str,
            'FT': str,
            'DPI': int,
            'PL': int,
            'PS': int,
            'SB': float
        }
        SVIM = {
            'ID': str,
            'SVTYPE': str,
            'END': int,
            'SVLEN': int,
            'SUPPORT': int,
            'STD_SPAN': float,
            'STD_POS': float,
            'STD_POS1': float,
            'STD_POS2': float,
            'ZMWS': int,
            'SEQS': str,
            'READS': str,
            'CUTPASTE': bool,
            'GT': bool,
            'DP': int,
            'AD': str,
            'CN': int
        }
        SVISIONPRO = {
            'ID': str,
            'END': int,
            'IT': str,
            'SVLEN': int,
            'BKPS': str,
            'BKPSIT': str,
            'SVTYPE': str,
            'SUPPORT': int,
            'VAF': float,
            'RNAMES': str,
            'GT': str,
            'DR': int,
            'DV': int,
            'PRECISE': bool
        }


class VariantCallTags:
    PASSED = 'passed'
    FAILED_FILTER = 'failed_filter'
    HOMOPOLYMER_REGION = 'homopolymer_region'
    NEARBY_EXCLUDED_VARIANT = 'nearby_excluded_variant'
    NEARBY_EXCLUDED_REGION = 'nearby_excluded_region'


class VariantFilterOperators:
    LESS_THAN = '<'
    LESS_THAN_OR_EQUAL_TO = '<='
    GREATER_THAN = '>'
    GREATER_THAN_OR_EQUAL_TO = '>='
    EQUALS = '=='
    NOT_EQUALS = '!='
    IN = 'in'


class VariantFilterQuantifiers:
    ALL = 'all'
    ANY = 'any'
    MEDIAN = 'median'
    AVERAGE = 'average'
    MIN = 'min'
    MAX = 'max'


class VariantFilterSampleTypes:
    CASE = 'case'
    CONTROL = 'control'


class VariantTypes:
    SINGLE_NUCLEOTIDE_VARIANT = 'SNV'
    MULTI_NUCLEOTIDE_VARIANT = 'MNV'
    INSERTION = 'INS'
    DELETION = 'DEL'
    INVERSION = 'INV'
    DUPLICATION = 'DUP'
    TRANSLOCATION = 'TRA'
    BREAKPOINT = 'BND'
    REFERENCE = 'REF'

    ALL = [
        SINGLE_NUCLEOTIDE_VARIANT,
        MULTI_NUCLEOTIDE_VARIANT,
        INSERTION,
        DELETION,
        INVERSION,
        DUPLICATION,
        TRANSLOCATION,
        BREAKPOINT
    ]

    FULL_NAMES = {
        REFERENCE: 'Reference',
        SINGLE_NUCLEOTIDE_VARIANT: 'Single-nucleotide Variant',
        INSERTION: 'Insertion',
        DELETION: 'Deletion',
        DUPLICATION: 'Duplication',
        INVERSION: 'Inversion',
        TRANSLOCATION: 'Translocation'
    }

    QueryTypeDictionary = {
        SINGLE_NUCLEOTIDE_VARIANT: [SINGLE_NUCLEOTIDE_VARIANT],
        MULTI_NUCLEOTIDE_VARIANT: [MULTI_NUCLEOTIDE_VARIANT],
        INSERTION: [DUPLICATION, INSERTION],
        DELETION: [DELETION],
        INVERSION: [BREAKPOINT, INVERSION, TRANSLOCATION],
        DUPLICATION: [DUPLICATION, INSERTION],
        TRANSLOCATION: [BREAKPOINT, INVERSION, TRANSLOCATION],
        BREAKPOINT: [BREAKPOINT, INVERSION, TRANSLOCATION],
    }

    class DuplicationSubtypes:
        TANDEM_DUPLICATION = 'DUP_TANDEM'
        SEGMENTAL_DUPLICATION = 'DUP_SEG'
        INTERSPERSED_DUPLICATION = 'DUP_INTERSPERSED'
        ALL = [
            TANDEM_DUPLICATION,
            SEGMENTAL_DUPLICATION,
            INTERSPERSED_DUPLICATION
        ]
