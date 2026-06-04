from vstolib.vcf import *
from ..data import get_data_path


def test_vcf2tsv_manta_1():
    vcf_file = get_data_path(name='vcf/manta-germline/colo829bl-dna-illumina_manta_outputs/candidateSmallIndels.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_manta_vcf(
        df_vcf=df_vcf,
        source_id='colo829bl',
        platform='illumina'
    )
    assert len(df_variants) <= 107907


def test_vcf2tsv_manta_2():
    vcf_file = get_data_path(name='vcf/manta-germline/colo829bl-dna-illumina_manta_outputs/candidateSV.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_manta_vcf(
        df_vcf=df_vcf,
        source_id='colo829bl',
        platform='illumina'
    )
    assert len(df_variants) <= 127162


def test_vcf2tsv_manta_3():
    vcf_file = get_data_path(name='vcf/manta-germline/hcc1395bl-dna-illumina_manta_outputs/candidateSmallIndels.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_manta_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395bl',
        platform='illumina'
    )
    assert len(df_variants) <= 119840


def test_vcf2tsv_manta_4():
    vcf_file = get_data_path(name='vcf/manta-germline/hcc1395bl-dna-illumina_manta_outputs/candidateSV.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_manta_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395bl',
        platform='illumina'
    )
    assert len(df_variants) <= 147495


def test_vcf2tsv_manta_5():
    vcf_file = get_data_path(name='vcf/manta-somatic/colo829-dna-illumina_manta_outputs/somaticSV.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_manta_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='illumina'
    )
    assert len(df_variants) >= 593


def test_vcf2tsv_manta_6():
    vcf_file = get_data_path(name='vcf/manta-somatic/hcc1395-dna-illumina_manta_outputs/somaticSV.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_manta_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='illumina'
    )
    assert len(df_variants) >= 1459
