from vstolib.constants import VariantFilterQuantifiers, VariantFilterOperators
from vstolib.main import filter, filter_excluded_regions, filter_homopolymeric_variants, score
from vstolib.variant_filter import VariantFilter
from .data import get_data_path


def test_filter_cutesv_variants_list(
        cutesv_variants_list,
        hg38_excluded_regions_list):
    variant_filters = []
    variant_filter_1 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='alternate_allele_read_count',
        operator=VariantFilterOperators.GREATER_THAN_OR_EQUAL_TO,
        value=2,
        sample_ids=['HG002']
    )
    variant_filter_2 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='filter',
        operator=VariantFilterOperators.EQUALS,
        value='PASS',
        sample_ids=['HG002']
    )
    variant_filters.append(variant_filter_1)
    variant_filters.append(variant_filter_2)
    variants_list_passed, variants_list_rejected = filter(
        variants_list=cutesv_variants_list,
        variant_filters=variant_filters,
        num_threads=1
    )
    variants_list_passed, variants_list_rejected = filter_excluded_regions(
        variants_list=variants_list_passed,
        excluded_regions_list=hg38_excluded_regions_list,
        excluded_regions_padding=1000,
        num_threads=1
    )


def test_filter_pbsv_variants_list(
        pbsv_variants_list,
        hg38_germline_variants_list,
        hg38_excluded_regions_list):
    variant_filters = []
    variant_filter_1 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='alternate_allele_read_count',
        operator=VariantFilterOperators.GREATER_THAN_OR_EQUAL_TO,
        value=2,
        sample_ids=['HG002']
    )
    variant_filter_2 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='filter',
        operator=VariantFilterOperators.EQUALS,
        value='PASS',
        sample_ids=['HG002']
    )
    variant_filters.append(variant_filter_1)
    variant_filters.append(variant_filter_2)
    variants_list_passed, variants_list_rejected = filter(
        variants_list=pbsv_variants_list,
        variant_filters=variant_filters,
        num_threads=1
    )
    variants_list_passed, variants_list_rejected = filter_excluded_regions(
        variants_list=variants_list_passed,
        excluded_regions_list=hg38_excluded_regions_list,
        excluded_regions_padding=1000,
        num_threads=1
    )


def test_filter_sniffles2_variants_list(
        sniffles2_variants_list,
        hg38_germline_variants_list,
        hg38_excluded_regions_list):
    variant_filter_1 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='alternate_allele_read_count',
        operator=VariantFilterOperators.GREATER_THAN_OR_EQUAL_TO,
        value=2,
        sample_ids=['HG002']
    )
    variant_filter_2 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='filter',
        operator=VariantFilterOperators.EQUALS,
        value='PASS',
        sample_ids=['HG002']
    )
    variant_filter_3 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='chromosome_1',
        operator=VariantFilterOperators.IN,
        value=['chr1'],
        sample_ids=['HG002']
    )
    variant_filter_4 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='chromosome_2',
        operator=VariantFilterOperators.IN,
        value=['chr1'],
        sample_ids=['HG002']
    )
    variant_filters = [variant_filter_1,
                       variant_filter_2,
                       variant_filter_3,
                       variant_filter_4]
    variants_list_passed, variants_list_rejected = filter_excluded_regions(
        variants_list=sniffles2_variants_list,
        excluded_regions_list=hg38_excluded_regions_list,
        excluded_regions_padding=1000,
        num_threads=1
    )


def test_filter_svim_variants_list(
        svim_variants_list,
        hg38_germline_variants_list,
        hg38_excluded_regions_list):
    variant_filters = []
    variant_filter_1 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='alternate_allele_read_count',
        operator=VariantFilterOperators.GREATER_THAN_OR_EQUAL_TO,
        value=2,
        sample_ids=['HG002']
    )
    variant_filter_2 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='filter',
        operator=VariantFilterOperators.EQUALS,
        value='PASS',
        sample_ids=['HG002']
    )
    variant_filters.append(variant_filter_1)
    variant_filters.append(variant_filter_2)
    variants_list_passed, variants_list_rejected = filter_excluded_regions(
        variants_list=svim_variants_list,
        excluded_regions_list=hg38_excluded_regions_list,
        excluded_regions_padding=1000,
        num_threads=1
    )


def test_filter_deepvariant_variants_list(
        deepvariant_variants_list,
        hg38_germline_variants_list,
        hg38_excluded_regions_list):
    variant_filters = []
    variant_filter_1 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='alternate_allele_read_count',
        operator=VariantFilterOperators.GREATER_THAN_OR_EQUAL_TO,
        value=2,
        sample_ids=['HG002']
    )
    variant_filter_2 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='filter',
        operator=VariantFilterOperators.EQUALS,
        value='PASS',
        sample_ids=['HG002']
    )
    variant_filters.append(variant_filter_1)
    variant_filters.append(variant_filter_2)
    variants_list_passed, variants_list_rejected = filter_excluded_regions(
        variants_list=deepvariant_variants_list,
        excluded_regions_list=hg38_excluded_regions_list,
        excluded_regions_padding=1000,
        num_threads=1
    )


def test_filter_gatk4_mutect2_variants_list(
        gatk4_mutect2_variants_list,
        hg38_germline_variants_list,
        hg38_excluded_regions_list):
    variant_filters = []
    variant_filter_1 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='alternate_allele_read_count',
        operator=VariantFilterOperators.GREATER_THAN_OR_EQUAL_TO,
        value=2,
        sample_ids=['hg002']
    )
    variant_filter_2 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='filter',
        operator=VariantFilterOperators.EQUALS,
        value='PASS',
        sample_ids=['hg002']
    )
    variant_filters.append(variant_filter_1)
    variant_filters.append(variant_filter_2)
    variants_list_passed, variants_list_rejected = filter_excluded_regions(
        variants_list=gatk4_mutect2_variants_list,
        excluded_regions_list=hg38_excluded_regions_list,
        excluded_regions_padding=1000,
        num_threads=1
    )


def test_filter_severus_variants_list(severus_variants_list):
    bam_file = get_data_path(name='hg38_tp53_tumor_long_read_dna.bam')
    variants_list = score(
        variants_list=severus_variants_list,
        bam_file=bam_file,
        window=1000,
        num_threads=1
    )
    variant_filter_1 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='position_1_average_alignment_score',
        operator=VariantFilterOperators.GREATER_THAN_OR_EQUAL_TO,
        value=30,
        sample_ids=['tumor']
    )
    variant_filters = [variant_filter_1]
    variants_list_passed, variants_list_rejected = filter(
        variants_list=variants_list,
        variant_filters=variant_filters,
        num_threads=1
    )