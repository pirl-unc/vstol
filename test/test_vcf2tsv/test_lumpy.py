from vstolib.vcf import *
from ..data import get_data_path


def test_vcf2tsv_lumpy_1():
    vcf_file = get_data_path(name='vcf/lumpy-germline/colo829bl-dna-illumina_lumpy.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_lumpy_vcf(
        df_vcf=df_vcf,
        source_id='colo829bl',
        platform='illumina'
    )
    assert len(df_variants) >= 13996


def test_vcf2tsv_lumpy_2():
    vcf_file = get_data_path(name='vcf/lumpy-germline/hcc1395bl-dna-illumina_lumpy.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_lumpy_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395bl',
        platform='illumina'
    )
    assert len(df_variants) >= 19724


def test_vcf2tsv_lumpy_3():
    vcf_file = get_data_path(name='vcf/lumpy-somatic/colo829-dna-illumina_lumpy.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_lumpy_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='illumina'
    )
    assert len(df_variants) >= 37444


def test_vcf2tsv_lumpy_4():
    vcf_file = get_data_path(name='vcf/lumpy-somatic/hcc1395-dna-illumina_lumpy.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_lumpy_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='illumina'
    )
    assert len(df_variants) >= 27444