from vstolib.main import overlap


def test_overlap_sniffles2_variants_list(
        sniffles2_variants_list,
        hg38_excluded_regions_list):
    variants_list = overlap(
        variants_list=sniffles2_variants_list,
        genomic_ranges_list=hg38_excluded_regions_list,
        num_threads=1
    )
