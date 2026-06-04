from vstolib.vcf import *
from ..data import get_data_path


def test_vcf2tsv_severus_1():
    vcf_file = get_data_path(name='vcf/severus/colo829-dna-pacbio_severus_outputs/somatic_SVs/severus_somatic.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_severus_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='pacbio'
    )
    assert len(df_variants) <= 146


def test_vcf2tsv_severus_2():
    vcf_file = get_data_path(name='vcf/severus/hcc1395-dna-pacbio_severus_outputs/somatic_SVs/severus_somatic.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_severus_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='pacbio'
    )
    assert len(df_variants) <= 2249