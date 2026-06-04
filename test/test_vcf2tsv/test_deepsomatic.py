from vstolib.vcf import *
from ..data import get_data_path


def test_vcf2tsv_deepsomatic_1():
    vcf_file = get_data_path(name='vcf/deepsomatic/colo829-dna-pacbio_deepsomatic.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_deepsomatic_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='pacbio'
    )
    assert len(df_variants) >= 1070231


def test_vcf2tsv_deepsomatic_2():
    vcf_file = get_data_path(name='vcf/deepsomatic/hcc1395-dna-pacbio_deepsomatic.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_deepsomatic_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='pacbio'
    )
    assert len(df_variants) >= 624539
