from vstolib.vcf import *
from ..data import get_data_path


def test_vcf2tsv_savana_1():
    vcf_file = get_data_path(name='vcf/savana/colo829-dna-pacbio_savana_classify_output.germline.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_savana_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='pacbio'
    )
    assert len(df_variants) == 37084


def test_vcf2tsv_savana_2():
    vcf_file = get_data_path(name='vcf/savana/colo829-dna-pacbio_savana_classify_output.somatic.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_savana_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='pacbio'
    )
    assert len(df_variants) == 554


def test_vcf2tsv_savana_3():
    vcf_file = get_data_path(name='vcf/savana/hcc1395-dna-pacbio_savana_classify_output.germline.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_savana_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='pacbio'
    )
    assert len(df_variants) == 31273


def test_vcf2tsv_savana_4():
    vcf_file = get_data_path(name='vcf/savana/hcc1395-dna-pacbio_savana_classify_output.somatic.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_savana_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='pacbio'
    )
    assert len(df_variants) == 3065