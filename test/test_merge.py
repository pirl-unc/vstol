import pandas as pd
from vstolib.main import merge


def test_merge_1():
    data = {
        'variant_id': [1, 2],
        'chromosome_1': ['chr17', 'chr17'],
        'position_1': [7673800, 7679800],
        'operation_1': ['D', 'D'],
        'chromosome_2': ['chr17', 'chr17'],
        'position_2': [7673801, 7679801],
        'operation_2': ['U', 'U'],
        'variant_size': [1000, 1000]
    }

    df_variants_1 = pd.DataFrame(data)

    data = {
        'variant_id': [100],
        'chromosome_1': ['chr17'],
        'position_1': [7673790],
        'operation_1': ['D'],
        'chromosome_2': ['chr17'],
        'position_2': [7673791],
        'operation_2': ['U'],
        'variant_size': [1000]
    }

    df_variants_2 = pd.DataFrame(data)

    df_intersecting = merge(
        df_variants_list=[df_variants_1, df_variants_2],
        match_both_positions=True,
        max_breakpoint_distance=1000,
        max_small_var_size=30,
        match_operation_types=True,
        num_processes=1
    )

    assert len(df_intersecting['group_id'].unique()) == 2


def test_merge_input_with_preexisting_group_id():
    # Regression: an input that already carries a `group_id` column (e.g. a
    # previously-merged or collapsed TSV) used to crash with
    # "All arrays must be of the same length" because `group_id` was appended
    # twice per row. The stale column must be dropped and fresh group_ids
    # assigned.
    df_variants_1 = pd.DataFrame({
        'variant_id': [1, 2],
        'chromosome_1': ['chr17', 'chr17'],
        'position_1': [7673800, 7679800],
        'operation_1': ['D', 'D'],
        'chromosome_2': ['chr17', 'chr17'],
        'position_2': [7673801, 7679801],
        'operation_2': ['U', 'U'],
        'variant_size': [1000, 1000],
        'group_id': [7, 7],  # stale group_id from an upstream merge/collapse
    })

    df_variants_2 = pd.DataFrame({
        'variant_id': [100],
        'chromosome_1': ['chr17'],
        'position_1': [7673790],
        'operation_1': ['D'],
        'chromosome_2': ['chr17'],
        'position_2': [7673791],
        'operation_2': ['U'],
        'variant_size': [1000],
    })

    df_merged = merge(
        df_variants_list=[df_variants_1, df_variants_2],
        match_both_positions=True,
        max_breakpoint_distance=1000,
        max_small_var_size=30,
        match_operation_types=True,
        num_processes=1
    )

    # Did not raise, and group_ids are freshly assigned (not the stale 7).
    assert len(df_merged) == 3
    assert df_merged['group_id'].nunique() == 2
    assert 7 not in set(df_merged['group_id'])