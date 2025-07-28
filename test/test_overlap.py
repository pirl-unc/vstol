import pandas as pd
from vstolib.main import overlap
from .data import get_data_path


def test_overlap_1():
    data = {
        'id': [1, 2],
        'chromosome_1': ['chr17', 'chr17'],
        'position_1': [7673800, 7679800],
        'operation_1': ['D', 'D'],
        'chromosome_2': ['chr17', 'chr17'],
        'position_2': [7673801, 7679801],
        'operation_2': ['U', 'U']
    }

    df_variants = pd.DataFrame(data)

    tsv_file = get_data_path(name='tsv/hg38_ucsc_gap_table.tsv')
    df_genomic_ranges = pd.read_csv(tsv_file, sep='\t')

    df_overlapping = overlap(
        df_variants=df_variants,
        df_genomic_ranges=df_genomic_ranges,
        buffer=0
    )

    assert len(df_overlapping) == 0
