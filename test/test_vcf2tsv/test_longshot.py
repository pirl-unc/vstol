from vstolib.vcf import *
from ..data import get_data_path


def test_vcf2tsv_longshot_1():
    vcf_file = get_data_path(name='vcf/longshot/colo829bl-dna-pacbio_longshot.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_longshot_vcf(
        df_vcf=df_vcf,
        source_id='colo829bl',
        platform='pacbio'
    )
    assert len(df_variants) >= 313386


def test_vcf2tsv_longshot_2():
    vcf_file = get_data_path(name='vcf/longshot/hcc1395bl-dna-pacbio_longshot.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_longshot_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395bl',
        platform='pacbio'
    )
    assert len(df_variants) >= 309800