from vstolib.vcf import *
from ..data import get_data_path


def test_vcf2tsv_sniffles2_1():
    vcf_file = get_data_path(name='vcf/sniffles2/colo829bl-dna-pacbio_sniffles2.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_sniffles2_vcf(
        df_vcf=df_vcf,
        source_id='colo829bl',
        platform='pacbio'
    )
    assert len(df_variants) >= 45136


def test_vcf2tsv_sniffles2_2():
    vcf_file = get_data_path(name='vcf/sniffles2/hcc1395bl-dna-pacbio_sniffles2.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_sniffles2_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395bl',
        platform='pacbio'
    )
    assert len(df_variants) >= 41397
