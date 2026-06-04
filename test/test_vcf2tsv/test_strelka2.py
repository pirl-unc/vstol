from vstolib.vcf import *
from ..data import get_data_path


def test_vcf2tsv_strelka2_1():
    vcf_file = get_data_path(name='vcf/strelka2-germline/colo829bl-dna-illumina_strelka2.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_strelka2_germline_vcf(
        df_vcf=df_vcf,
        source_id='colo829bl',
        platform='illumina'
    )
    assert len(df_variants) >= 383990


def test_vcf2tsv_strelka2_2():
    vcf_file = get_data_path(name='vcf/strelka2-germline/hcc1395bl-dna-illumina_strelka2.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_strelka2_germline_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395bl',
        platform='illumina'
    )
    assert len(df_variants) >= 377847


def test_vcf2tsv_strelka2_3():
    vcf_file = get_data_path(name='vcf/strelka2-somatic/colo829-dna-illumina_strelka2_indels.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_strelka2_somatic_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='illumina'
    )
    assert len(df_variants) >= 31790


def test_vcf2tsv_strelka2_4():
    vcf_file = get_data_path(name='vcf/strelka2-somatic/colo829-dna-illumina_strelka2_snvs.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_strelka2_somatic_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='illumina'
    )
    assert len(df_variants) >= 169575


def test_vcf2tsv_strelka2_5():
    vcf_file = get_data_path(name='vcf/strelka2-somatic/hcc1395-dna-illumina_strelka2_indels.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_strelka2_somatic_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='illumina'
    )
    assert len(df_variants) >= 19655


def test_vcf2tsv_strelka2_6():
    vcf_file = get_data_path(name='vcf/strelka2-somatic/hcc1395-dna-illumina_strelka2_snvs.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_strelka2_somatic_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='illumina'
    )
    assert len(df_variants) >= 169575
