import pytest
from .data import get_data_path
from vstolib.constants import *
from vstolib.genomic_ranges_list import GenomicRangesList
from vstolib.main import vcf2tsv
from vstolib.variants_list import VariantsList
from vstolib.vcf.common import read_vcf_file


@pytest.fixture
def cutesv_variants_list():
    vcf_file = get_data_path(name='hg002_cutesv.vcf')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    variants_list = vcf2tsv(
        df_vcf=df_vcf,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.CUTESV,
        sequencing_platform='pacbio'
    )
    return variants_list

@pytest.fixture
def deepvariant_variants_list():
    vcf_file = get_data_path(name='hg002_deepvariant.vcf')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    variants_list = vcf2tsv(
        df_vcf=df_vcf,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.DEEPVARIANT,
        sequencing_platform='pacbio'
    )
    return variants_list

@pytest.fixture
def delly2_somatic_variants_list():
    vcf_file = get_data_path(name='hg002_hg001_delly2.vcf')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    variants_list = vcf2tsv(
        df_vcf=df_vcf,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.DELLY2_SOMATIC,
        sequencing_platform='ilmn'
    )
    return variants_list

@pytest.fixture
def gatk4_mutect2_variants_list():
    vcf_file = get_data_path(name='hg002_hg001_gatk4_mutect2.vcf')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    variants_list = vcf2tsv(
        df_vcf=df_vcf,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.GATK4_MUTECT2,
        sequencing_platform='ilmn'
    )
    return variants_list

@pytest.fixture
def lumpy_somatic_variants_list():
    vcf_file = get_data_path(name='hg002_hg001_lumpy.vcf')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    variants_list = vcf2tsv(
        df_vcf=df_vcf,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.LUMPY_SOMATIC,
        sequencing_platform='pacbio'
    )
    return variants_list

@pytest.fixture
def pbsv_variants_list():
    vcf_file = get_data_path(name='hg002_pbsv.vcf')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    variants_list = vcf2tsv(
        df_vcf=df_vcf,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.PBSV,
        sequencing_platform='pacbio'
    )
    return variants_list

@pytest.fixture
def sniffles2_variants_list():
    vcf_file = get_data_path(name='hg002_sniffles2.vcf')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    variants_list = vcf2tsv(
        df_vcf=df_vcf,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.SNIFFLES2,
        sequencing_platform='pacbio'
    )
    return variants_list

@pytest.fixture
def strelka2_somatic_indels_variants_list():
    vcf_file = get_data_path(name='hg002_hg001_strelka2_indels.vcf')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    variants_list = vcf2tsv(
        df_vcf=df_vcf,
        source_id='hg002',
        case_id='hg002',
        control_id='hg001',
        variant_calling_method=VariantCallingMethods.STRELKA2_SOMATIC,
        sequencing_platform='pacbio'
    )
    return variants_list

@pytest.fixture
def strelka2_somatic_snvs_variants_list():
    vcf_file = get_data_path(name='hg002_hg001_strelka2_snvs.vcf')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    variants_list = vcf2tsv(
        df_vcf=df_vcf,
        source_id='hg002',
        case_id='hg002',
        control_id='hg001',
        variant_calling_method=VariantCallingMethods.STRELKA2_SOMATIC,
        sequencing_platform='pacbio'
    )
    return variants_list

@pytest.fixture
def svim_variants_list():
    vcf_file = get_data_path(name='hg002_svim.vcf')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    variants_list = vcf2tsv(
        df_vcf=df_vcf,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.SVIM,
        sequencing_platform='pacbio'
    )
    return variants_list

@pytest.fixture
def hg38_germline_variants_list():
    germline_variants_tsv_file = get_data_path(name='audano_et_al_cell_2019_sv_list.tsv')
    germline_variants_list = VariantsList.read_tsv_file(tsv_file=germline_variants_tsv_file)
    return germline_variants_list

@pytest.fixture
def hg38_excluded_regions_list():
    excluded_regions_tsv_file = get_data_path(name='hg38_ucsc_gap_table.tsv')
    excluded_regions_list = GenomicRangesList.read_tsv_file(tsv_file=excluded_regions_tsv_file)
    return excluded_regions_list

