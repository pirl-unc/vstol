from vstolib.vcf import *
from ..data import get_data_path


def test_vcf2tsv_cutesv_1():
    vcf_file = get_data_path(name='vcf/cutesv/colo829bl-dna-pacbio_cutesv.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_cutesv_vcf(
        df_vcf=df_vcf,
        source_id='colo829bl',
        platform='pacbio'
    )
    assert len(df_variants) >= 2598


def test_vcf2tsv_cutesv_2():
    vcf_file = get_data_path(name='vcf/cutesv/hcc1395bl-dna-pacbio_cutesv.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_cutesv_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395bl',
        platform='pacbio'
    )
    assert len(df_variants) >= 2575
