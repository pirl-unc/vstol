from .conftest import *
from vstolib.main import compare
from vstolib.variant import Variant
from vstolib.variant_call import VariantCall
from vstolib.variants_list import VariantsList


def test_compare_1():
    variant_call_1 = VariantCall(
        id='variant_call_1',
        sample_id='sample1',
        chromosome_1='chr1',
        position_1=100,
        chromosome_2='chr1',
        position_2=100,
        variant_type=VariantTypes.SINGLE_NUCLEOTIDE_VARIANT,
        reference_allele='C',
        alternate_allele='A'
    )
    variant_call_2 = VariantCall(
        id='variant_call_2',
        sample_id='sample1',
        chromosome_1='chr1',
        position_1=1000,
        chromosome_2='chr1',
        position_2=1000,
        variant_type=VariantTypes.SINGLE_NUCLEOTIDE_VARIANT,
        reference_allele='T',
        alternate_allele='A'
    )
    variant_call_3 = VariantCall(
        id='variant_call_3',
        sample_id='sample1',
        chromosome_1='chr1',
        position_1=10000,
        chromosome_2='chr1',
        position_2=10000,
        variant_type=VariantTypes.SINGLE_NUCLEOTIDE_VARIANT,
        reference_allele='C',
        alternate_allele='A'
    )

    variant_1 = Variant(id='variant_1')
    variant_2 = Variant(id='variant_2')
    variant_3 = Variant(id='variant_3')
    variant_1.add_variant_call(variant_call=variant_call_1)
    variant_2.add_variant_call(variant_call=variant_call_2)
    variant_3.add_variant_call(variant_call=variant_call_3)

    variants_list_1 = VariantsList()
    variants_list_2 = VariantsList()
    variants_list_1.add_variant(variant=variant_1)
    variants_list_1.add_variant(variant=variant_2)
    variants_list_2.add_variant(variant=variant_1)
    variants_list_2.add_variant(variant=variant_3)

    vl_shared, vl_a_only, vl_b_only = compare(
        a=variants_list_1,
        b=variants_list_2,
        num_threads=1,
        max_neighbor_distance=100,
        match_all_breakpoints=True,
        match_variant_types=True,
        min_ins_size_overlap=0.5,
        min_del_size_overlap=0.5
    )

    assert(len(vl_shared.variants) == 1)
    assert(len(vl_a_only.variants) == 1)
    assert(len(vl_b_only.variants) == 1)
    assert(vl_shared.variants[0].variant_calls[0].id == 'variant_call_1')
    assert(vl_a_only.variants[0].variant_calls[0].id == 'variant_call_2')
    assert(vl_b_only.variants[0].variant_calls[0].id == 'variant_call_3')

