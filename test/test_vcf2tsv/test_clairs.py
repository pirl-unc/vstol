from vstolib.vcf import *
from ..data import get_data_path


def test_vcf2tsv_clairs_snv_1():
    vcf_file = get_data_path(name='vcf/clairs/colo829-dna-pacbio_clairs_outputs/snv.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_clairs_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='pacbio'
    )
    assert len(df_variants) == 54832


def test_vcf2tsv_clairs_snv_2():
    vcf_file = get_data_path(name='vcf/clairs/hcc1395-dna-pacbio_clairs_outputs/snv.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_clairs_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='pacbio'
    )
    assert len(df_variants) == 139851


def test_vcf2tsv_clairs_indel_1():
    vcf_file = get_data_path(name='vcf/clairs/colo829-dna-pacbio_clairs_outputs/indel.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_clairs_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='pacbio'
    )
    assert len(df_variants) == 180697


def test_vcf2tsv_clairs_indel_2():
    vcf_file = get_data_path(name='vcf/clairs/hcc1395-dna-pacbio_clairs_outputs/indel.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_clairs_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='pacbio'
    )
    assert len(df_variants) == 83012
