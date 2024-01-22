from vstolib.main import diff


def test_diff_sniffles2_variants_list(
        sniffles2_variants_list,
        cutesv_variants_list,
        pbsv_variants_list,
        svim_variants_list):
    sniffles2_variants_list_diff = diff(
        target_variants_list=sniffles2_variants_list,
        query_variants_lists=[
            cutesv_variants_list,
            pbsv_variants_list,
            svim_variants_list
        ],
        num_threads=1,
        max_neighbor_distance=10
    )
