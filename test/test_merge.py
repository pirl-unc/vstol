from .conftest import *
from vstolib.main import merge
from vstolib.variant import Variant
from vstolib.variant_call import VariantCall
from vstolib.variants_list import VariantsList


def test_merge_1(cutesv_variants_list,
                deepvariant_variants_list,
                pbsv_variants_list,
                sniffles2_variants_list,
                svim_variants_list):
    variants_list = merge(
        variants_lists=[cutesv_variants_list,
                        deepvariant_variants_list,
                        pbsv_variants_list,
                        sniffles2_variants_list,
                        svim_variants_list],
        num_threads=1,
        max_neighbor_distance=10,
        match_all_breakpoints=True,
        match_variant_types=True
    )
    print(variants_list.size)


def test_merge_2():
    variant_call_11 = VariantCall(
        id='variant_call_11',
        sample_id='sample1',
        chromosome_1='chr1',
        position_1=100,
        chromosome_2='chr5',
        position_2=500,
        variant_type=VariantTypes.TRANSLOCATION,
        reference_allele='T',
        alternate_allele=''
    )
    variant_call_12 = VariantCall(
        id='variant_call_12',
        sample_id='sample1',
        chromosome_1='chr1',
        position_1=111,
        chromosome_2='chr1',
        position_2=120,
        variant_type=VariantTypes.DELETION,
        reference_allele='TACTGATCGA',
        alternate_allele=''
    )
    variant_call_13 = VariantCall(
        id='variant_call_13',
        sample_id='sample1',
        chromosome_1='chr7',
        position_1=700,
        chromosome_2='chr8',
        position_2=800,
        variant_type=VariantTypes.TRANSLOCATION,
        reference_allele='A',
        alternate_allele=''
    )
    variant_call_14 = VariantCall(
        id='variant_call_14',
        sample_id='sample1',
        chromosome_1='chr10',
        position_1=9500,
        chromosome_2='chr10',
        position_2=9600,
        variant_type=VariantTypes.DELETION,
        reference_allele='',
        alternate_allele=''
    )

    variant_call_21 = VariantCall(
        id='variant_call_21',
        sample_id='sample2',
        chromosome_1='chr1',
        position_1=101,
        chromosome_2='chr16',
        variant_type=VariantTypes.TRANSLOCATION,
        position_2=1000,
        reference_allele='C',
        alternate_allele=''
    )
    variant_call_22 = VariantCall(
        id='variant_call_22',
        sample_id='sample2',
        chromosome_1='chr10',
        position_1=9505,
        chromosome_2='chr10',
        position_2=9605,
        variant_type=VariantTypes.DELETION,
        reference_allele='',
        alternate_allele=''
    )

    variant_11 = Variant(id='variant_11')
    variant_12 = Variant(id='variant_12')
    variant_13 = Variant(id='variant_13')
    variant_14 = Variant(id='variant_14')
    variant_21 = Variant(id='variant_21')
    variant_22 = Variant(id='variant_22')
    variant_11.add_variant_call(variant_call=variant_call_11)
    variant_12.add_variant_call(variant_call=variant_call_12)
    variant_13.add_variant_call(variant_call=variant_call_13)
    variant_14.add_variant_call(variant_call=variant_call_14)
    variant_21.add_variant_call(variant_call=variant_call_21)
    variant_22.add_variant_call(variant_call=variant_call_22)

    variants_list_1 = VariantsList()
    variants_list_2 = VariantsList()
    variants_list_1.add_variant(variant=variant_11)
    variants_list_1.add_variant(variant=variant_12)
    variants_list_1.add_variant(variant=variant_13)
    variants_list_1.add_variant(variant=variant_14)
    variants_list_2.add_variant(variant=variant_21)
    variants_list_2.add_variant(variant=variant_22)

    variants_list = merge(
        variants_lists=[variants_list_1, variants_list_2],
        num_threads=1,
        max_neighbor_distance=1000,
        match_all_breakpoints=True,
        match_variant_types=True
    )
    for variant in variants_list.variants:
        if variant.num_variant_calls == 2:
            assert 'variant_call_14' in variant.variant_call_ids
            assert 'variant_call_22' in variant.variant_call_ids
