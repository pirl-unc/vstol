from vstolib.vcf import *
from ..data import get_data_path


def test_vcf2tsv_svaba_1():
    vcf_file = get_data_path(name='vcf/svaba/colo829-dna-illumina_svaba_outputs/colo829-dna-illumina.svaba.somatic.indel.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_svaba_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='illumina'
    )
    assert len(df_variants) >= 2214


def test_vcf2tsv_svaba_2():
    vcf_file = get_data_path(name='vcf/svaba/colo829-dna-illumina_svaba_outputs/colo829-dna-illumina.svaba.somatic.sv.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_svaba_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='illumina'
    )
    assert len(df_variants) >= 1886


def test_vcf2tsv_svaba_3():
    vcf_file = get_data_path(name='vcf/svaba/hcc1395-dna-illumina_svaba_outputs/hcc1395-dna-illumina.svaba.somatic.indel.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_svaba_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='illumina'
    )
    assert len(df_variants) >= 15807


def test_vcf2tsv_svaba_4():
    vcf_file = get_data_path(name='vcf/svaba/hcc1395-dna-illumina_svaba_outputs/hcc1395-dna-illumina.svaba.somatic.sv.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_svaba_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='illumina'
    )
    assert len(df_variants) >= 824

