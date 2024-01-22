from .conftest import *
from vstolib.main import intersect


def test_intersect(pbsv_variants_list,
                   sniffles2_variants_list):
    variants_list = intersect(
        variants_lists=[pbsv_variants_list,
                        sniffles2_variants_list],
        num_threads=1,
        max_neighbor_distance=10
    )
    print(variants_list.size)
