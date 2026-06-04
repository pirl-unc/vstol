from vstolib.vcf import *
from ..data import get_data_path


def test_vcf2tsv_mutect2_1():
    vcf_file = get_data_path(name='vcf/mutect2/colo829-dna-illumina_gatk4-mutect2.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_mutect2_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='illumina'
    )
    assert len(df_variants) >= 306914


def test_vcf2tsv_mutect2_2():
    vcf_file = get_data_path(name='vcf/mutect2/hcc1395-dna-illumina_gatk4-mutect2.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_mutect2_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='illumina'
    )
    assert len(df_variants) >= 256498