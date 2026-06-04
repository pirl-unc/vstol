from vstolib.vcf import *
from ..data import get_data_path


def test_vcf2tsv_dysgu_germline_1():
    vcf_file = get_data_path(name='vcf/dysgu-germline/colo829bl-dna-pacbio_dysgu_germline_svs.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_dysgu_vcf(
        df_vcf=df_vcf,
        source_id='colo829bl',
        platform='pacbio'
    )
    assert len(df_variants) >= 40682


def test_vcf2tsv_dysgu_germline_2():
    vcf_file = get_data_path(name='vcf/dysgu-germline/hcc1395bl-dna-pacbio_dysgu_germline_svs.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_dysgu_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395bl',
        platform='pacbio'
    )
    assert len(df_variants) >= 38568


def test_vcf2tsv_dysgu_germline_3():
    vcf_file = get_data_path(name='vcf/dysgu-germline/colo829bl-dna-illumina_dysgu_germline_svs.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_dysgu_vcf(
        df_vcf=df_vcf,
        source_id='colo829bl',
        platform='illumina'
    )
    assert len(df_variants) >= 22655


def test_vcf2tsv_dysgu_germline_4():
    vcf_file = get_data_path(name='vcf/dysgu-germline/hcc1395bl-dna-illumina_dysgu_germline_svs.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_dysgu_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395bl',
        platform='illumina'
    )
    assert len(df_variants) >= 25951


def test_vcf2tsv_dysgu_somatic_1():
    vcf_file = get_data_path(name='vcf/dysgu-somatic/colo829-dna-pacbio_dysgu_somatic_svs.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_dysgu_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='pacbio'
    )
    assert len(df_variants) >= 11


def test_vcf2tsv_dysgu_somatic_2():
    vcf_file = get_data_path(name='vcf/dysgu-somatic/hcc1395-dna-pacbio_dysgu_somatic_svs.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_dysgu_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='pacbio'
    )
    assert len(df_variants) >= 103


def test_vcf2tsv_dysgu_somatic_3():
    vcf_file = get_data_path(name='vcf/dysgu-somatic/colo829-dna-illumina_dysgu_somatic_svs.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_dysgu_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='illumina'
    )
    assert len(df_variants) >= 1990


def test_vcf2tsv_dysgu_somatic_4():
    vcf_file = get_data_path(name='vcf/dysgu-somatic/hcc1395-dna-illumina_dysgu_somatic_svs.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_dysgu_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='illumina'
    )
    assert len(df_variants) >= 1999