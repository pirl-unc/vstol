from .conftest import *


def test_vcf2tsv_clairs(clairs_variants_list):
    print(clairs_variants_list.size)


def test_vcf2tsv_cutesv(cutesv_variants_list):
    print(cutesv_variants_list.size)


def test_vcf2tsv_deepvariant(deepvariant_variants_list):
    print(deepvariant_variants_list.size)


def test_vcf2tsv_delly2_somatic(delly2_somatic_variants_list):
    print(delly2_somatic_variants_list.size)


def test_vcf2tsv_gatk4_mutect2(gatk4_mutect2_variants_list):
    print(gatk4_mutect2_variants_list.size)


def test_vcf2tsv_lumpy_somatic(lumpy_somatic_variants_list):
    print(lumpy_somatic_variants_list.size)


def test_vcf2tsv_manta_somatic(manta_somatic_variants_list):
    print(manta_somatic_variants_list.size)


def test_vcf2tsv_pbsv(pbsv_variants_list):
    print(pbsv_variants_list.size)


def test_vcf2tsv_savana(savana_variants_list):
    print(savana_variants_list.size)


def test_vcf2tsv_severus(severus_variants_list):
    print(severus_variants_list.size)


def test_vcf2tsv_sniffles2(sniffles2_variants_list):
    print(sniffles2_variants_list.size)


def test_vcf2tsv_strelka2_somatic_indels(strelka2_somatic_indels_variants_list):
    print(strelka2_somatic_indels_variants_list.size)


def test_vcf2tsv_strelka2_somatic_snvs(strelka2_somatic_snvs_variants_list):
    print(strelka2_somatic_snvs_variants_list.size)


def test_vcf2tsv_svim(svim_variants_list):
    print(svim_variants_list.size)


def test_vcf2tsv_svisionpro(svisionpro_variants_list):
    print(svisionpro_variants_list.size)
