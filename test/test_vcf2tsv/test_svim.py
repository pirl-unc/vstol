from vstolib.vcf import *
from ..data import get_data_path


def test_vcf2tsv_svim_1():
    vcf_file = get_data_path(name='vcf/svim/colo829bl-dna-pacbio_svim.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_svim_vcf(
        df_vcf=df_vcf,
        source_id='colo829bl',
        platform='pacbio'
    )
    assert len(df_variants) >= 140359


def test_vcf2tsv_svim_2():
    vcf_file = get_data_path(name='vcf/svim/hcc1395bl-dna-pacbio_svim.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_svim_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395bl',
        platform='pacbio'
    )
    assert len(df_variants) >= 101328
