from vstolib.vcf import *
from ..data import get_data_path


def test_vcf2tsv_haplotypecaller_1():
    vcf_file = get_data_path(name='vcf/haplotypecaller/colo829bl-dna-illumina_gatk4-haplotypecaller_chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_haplotypecaller_vcf(
        df_vcf=df_vcf,
        source_id='colo829bl',
        platform='pacbio'
    )
    assert len(df_variants) >= 366598


def test_vcf2tsv_haplotypecaller_2():
    vcf_file = get_data_path(name='vcf/haplotypecaller/hcc1395bl-dna-illumina_gatk4-haplotypecaller_chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_haplotypecaller_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395bl',
        platform='pacbio'
    )
    assert len(df_variants) >= 366272
