from vstolib.vcf import *
from ..data import get_data_path


def test_vcf2tsv_delly2_1():
    vcf_file = get_data_path(name='vcf/delly2-germline/colo829bl-dna-illumina_delly2.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_delly2_vcf(
        df_vcf=df_vcf,
        source_id='colo829bl',
        platform='illumina'
    )
    assert len(df_variants) >= 37904


def test_vcf2tsv_delly2_2():
    vcf_file = get_data_path(name='vcf/delly2-germline/hcc1395bl-dna-illumina_delly2.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_delly2_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395bl',
        platform='illumina'
    )
    assert len(df_variants) >= 48786


def test_vcf2tsv_delly2_3():
    vcf_file = get_data_path(name='vcf/delly2-somatic/colo829-dna-illumina_delly2_filtered.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_delly2_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='illumina'
    )
    assert len(df_variants) >= 269


def test_vcf2tsv_delly2_4():
    vcf_file = get_data_path(name='vcf/delly2-somatic/hcc1395-dna-illumina_delly2_filtered.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_delly2_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='illumina'
    )
    assert len(df_variants) >= 1410

