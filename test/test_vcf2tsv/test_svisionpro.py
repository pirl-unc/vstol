from vstolib.vcf import *
from ..data import get_data_path


def test_vcf2tsv_svisionpro_1():
    vcf_file = get_data_path(name='vcf/svisionpro/colo829-dna-pacbio.svision_pro_v1.8.s3.somatic_s3.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_svisionpro_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='pacbio'
    )
    assert len(df_variants) >= 215


def test_vcf2tsv_svisionpro_2():
    vcf_file = get_data_path(name='vcf/svisionpro/hcc1395-dna-pacbio.svision_pro_v1.8.s3.somatic_s3.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_svisionpro_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='pacbio'
    )
    assert len(df_variants) >= 2102
