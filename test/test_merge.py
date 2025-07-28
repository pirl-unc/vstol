import pandas as pd
from vstolib.main import merge


def test_merge_1():
    data = {
        'id': [1, 2],
        'chromosome_1': ['chr17', 'chr17'],
        'position_1': [7673800, 7679800],
        'operation_1': ['D', 'D'],
        'chromosome_2': ['chr17', 'chr17'],
        'position_2': [7673801, 7679801],
        'operation_2': ['U', 'U']
    }

    df_variants_1 = pd.DataFrame(data)

    data = {
        'id': [100],
        'chromosome_1': ['chr17'],
        'position_1': [7673790],
        'operation_1': ['D'],
        'chromosome_2': ['chr17'],
        'position_2': [7673791],
        'operation_2': ['U']
    }

    df_variants_2 = pd.DataFrame(data)

    df_intersecting = merge(
        df_variants_list=[df_variants_1, df_variants_2],
        match_both_positions=True,
        max_breakpoint_distance=1000,
        match_operation_types=True,
        num_processes=1
    )

    assert len(df_intersecting['group_id'].unique()) == 2