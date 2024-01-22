from .conftest import *
from vstolib.main import merge


def test_merge(cutesv_variants_list,
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
        max_neighbor_distance=10
    )
    print(variants_list.size)
