from vstolib.vcf import *
from ..data import get_data_path


def test_vcf2tsv_nanocaller_1():
    vcf_file = get_data_path(name='vcf/nanocaller/colo829bl-dna-pacbio_nanocaller_outputs/colo829bl-dna-pacbio_nanocaller.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_nanocaller_vcf(
        df_vcf=df_vcf,
        source_id='colo829bl',
        platform='pacbio'
    )
    assert len(df_variants) <= 2093423


def test_vcf2tsv_nanocaller_2():
    vcf_file = get_data_path(name='vcf/nanocaller/hcc1395bl-dna-pacbio_nanocaller_outputs/hcc1395bl-dna-pacbio_nanocaller.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_nanocaller_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395bl',
        platform='pacbio'
    )
    assert len(df_variants) <= 2061733