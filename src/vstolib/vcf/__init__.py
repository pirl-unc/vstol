from .common import read_vcf_file

from .clairs import parse_clairs_vcf
from .cutesv import parse_cutesv_vcf
from .dbsnp import parse_dbsnp_vcf
from .deepvariant import parse_deepvariant_vcf
from .delly2 import parse_delly2_somatic_vcf
from .gatk4_mutect2 import parse_gatk4_mutect2_vcf
from .lumpy import parse_lumpy_somatic_vcf
from .manta import parse_manta_somatic_vcf
from .pbsv import parse_pbsv_vcf
from .savana import parse_savana_vcf
from .severus import parse_severus_vcf
from .sniffles2 import parse_sniffles2_vcf
from .strelka2 import parse_strelka2_vcf
from .svim import parse_svim_vcf
from .svisionpro import parse_svisionpro_vcf