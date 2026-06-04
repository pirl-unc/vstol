from vstolib.vcf import *
from ..data import get_data_path


def test_vcf2tsv_nanomonsv_1():
    vcf_file = get_data_path(name='vcf/nanomonsv/colo829-dna-pacbio_nanomonsv_outputs/colo829-dna-pacbio.nanomonsv.result.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_nanomonsv_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='pacbio'
    )
    assert len(df_variants) >= 160


def test_vcf2tsv_nanomonsv_2():
    vcf_file = get_data_path(name='vcf/nanomonsv/hcc1395-dna-pacbio_nanomonsv_outputs/hcc1395-dna-pacbio.nanomonsv.result.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_nanomonsv_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='pacbio'
    )
    assert len(df_variants) >= 2401