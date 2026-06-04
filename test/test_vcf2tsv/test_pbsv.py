from vstolib.vcf import *
from ..data import get_data_path


def test_vcf2tsv_pbsv_1():
    vcf_file = get_data_path(name='vcf/pbsv/colo829bl-dna-pacbio_pbsv.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_pbsv_vcf(
        df_vcf=df_vcf,
        source_id='colo829bl',
        platform='pacbio'
    )
    assert len(df_variants) >= 53135


def test_vcf2tsv_pbsv_2():
    vcf_file = get_data_path(name='vcf/pbsv/hcc1395bl-dna-pacbio_pbsv.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_pbsv_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395bl',
        platform='pacbio'
    )
    assert len(df_variants) >= 51547
