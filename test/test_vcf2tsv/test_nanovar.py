from vstolib.vcf import *
from ..data import get_data_path


def test_vcf2tsv_nanovar_1():
    vcf_file = get_data_path(name='vcf/nanovar/colo829bl-dna-pacbio_nanovar_outputs/colo829bl-dna-pacbio_minimap2_mdtagged_sorted.nanovar.pass.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_nanovar_vcf(
        df_vcf=df_vcf,
        source_id='colo829bl',
        platform='pacbio'
    )
    assert len(df_variants) >= 44233


def test_vcf2tsv_nanovar_2():
    vcf_file = get_data_path(name='vcf/nanovar/hcc1395bl-dna-pacbio_nanovar_outputs/hcc1395bl-dna-pacbio_minimap2_mdtagged_sorted.nanovar.pass.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_nanovar_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395bl',
        platform='pacbio'
    )
    assert len(df_variants) >= 41932