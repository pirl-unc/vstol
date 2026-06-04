from vstolib.vcf import *
from ..data import get_data_path


def test_vcf2tsv_octopus_1():
    vcf_file = get_data_path(name='vcf/octopus-germline/colo829bl-dna-illumina_octopus.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_octopus_vcf(
        df_vcf=df_vcf,
        source_id='colo829bl',
        platform='pacbio'
    )
    assert len(df_variants) >= 383800


def test_vcf2tsv_octopus_2():
    vcf_file = get_data_path(name='vcf/octopus-germline/hcc1395bl-dna-illumina_octopus.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_octopus_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395bl',
        platform='pacbio'
    )
    assert len(df_variants) >= 372544


def test_vcf2tsv_octopus_3():
    vcf_file = get_data_path(name='vcf/octopus-somatic/hcc1395-dna-illumina_octopus.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_octopus_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='pacbio'
    )
    assert len(df_variants) >= 382174
