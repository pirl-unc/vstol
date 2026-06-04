from vstolib.vcf import *
from ..data import get_data_path


def test_vcf2tsv_clair3_1():
    vcf_file = get_data_path(name='vcf/clair3/colo829bl-dna-pacbio_clair3_outputs/merge_output.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_clair3_vcf(
        df_vcf=df_vcf,
        source_id='colo829bl',
        platform='pacbio'
    )
    assert len(df_variants) >= 393441


def test_vcf2tsv_clair3_2():
    vcf_file = get_data_path(name='vcf/clair3/hcc1395bl-dna-pacbio_clair3_outputs/merge_output.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_clair3_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395bl',
        platform='pacbio'
    )
    assert len(df_variants) >= 392278


def test_vcf2tsv_clair3_3():
    vcf_file = get_data_path(name='vcf/clair3/colo829bl-dna-illumina_clair3_outputs/merge_output.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_clair3_vcf(
        df_vcf=df_vcf,
        source_id='colo829bl',
        platform='illumina'
    )
    assert len(df_variants) >= 373647


def test_vcf2tsv_clair3_4():
    vcf_file = get_data_path(name='vcf/clair3/hcc1395bl-dna-illumina_clair3_outputs/merge_output.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_clair3_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395bl',
        platform='illumina'
    )
    assert len(df_variants) >= 364639


