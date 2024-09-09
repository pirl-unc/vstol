from vstolib.main import score
from .data import get_data_path


def test_score(severus_variants_list):
    bam_file = get_data_path(name='hg38_tp53_tumor_long_read_dna.bam')
    variants_list = score(
        variants_list=severus_variants_list,
        bam_file=bam_file,
        window=1000,
        num_threads=1
    )
