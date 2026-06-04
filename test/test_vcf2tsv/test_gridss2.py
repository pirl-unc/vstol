from vstolib.vcf import *

from src.vstolib.vcf import parse_gridss2_vcf
from ..data import get_data_path


def test_vcf2tsv_gridss2_1():
    vcf_file = get_data_path(name='vcf/gridss2-germline/colo829bl-dna-illumina_gridss.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_gridss2_vcf(
        df_vcf=df_vcf,
        source_id='colo829bl',
        platform='pacbio'
    )
    assert len(df_variants) <= 41340


def test_vcf2tsv_gridss2_2():
    vcf_file = get_data_path(name='vcf/gridss2-germline/hcc1395bl-dna-illumina_gridss.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_gridss2_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395bl',
        platform='pacbio'
    )
    assert len(df_variants) <= 48153


def test_vcf2tsv_gridss2_3():
    vcf_file = get_data_path(name='vcf/gridss2-somatic/colo829-dna-illumina_gridss_high_confidence_somatic.vcf.bgz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_gridss2_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='pacbio'
    )
    assert len(df_variants) >= 211


def test_vcf2tsv_gridss2_4():
    vcf_file = get_data_path(name='vcf/gridss2-somatic/hcc1395-dna-illumina_gridss_high_confidence_somatic.vcf.bgz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_gridss2_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='pacbio'
    )
    assert len(df_variants) >= 2404