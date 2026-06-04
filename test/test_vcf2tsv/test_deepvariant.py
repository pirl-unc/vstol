from vstolib.vcf import *
from ..data import get_data_path


def test_vcf2tsv_deepvariant_1():
    vcf_file = get_data_path(name='vcf/deepvariant/colo829bl-dna-pacbio_deepvariant.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_deepvariant_vcf(
        df_vcf=df_vcf,
        source_id='colo829bl',
        platform='pacbio'
    )
    assert len(df_variants) >= 797988


def test_vcf2tsv_deepvariant_2():
    vcf_file = get_data_path(name='vcf/deepvariant/hcc1395bl-dna-pacbio_deepvariant.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_deepvariant_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395bl',
        platform='pacbio'
    )
    assert len(df_variants) >= 787461
