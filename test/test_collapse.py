import pandas as pd
import pytest
from vstolib.main import collapse


def _row(
        group_id,
        sample_id,
        method='savana',
        position_1=1000,
        position_2=2000,
        p1_alt=10,
        p2_alt=10,
        quality=50.0,
        variant_size=1000,
        chromosome='chr1',
        variant_type='DEL',
        chromosome_1=None,
        chromosome_2=None,
        variant_id=None,
        sequence='',
        filter_value='PASS',
):
    """Build a single fully-populated merged-TSV row dict.

    `chromosome` sets both breakends by default; pass `chromosome_1` /
    `chromosome_2` explicitly for interchromosomal events (translocations).
    `variant_id` defaults to `group_id`; pass it explicitly to exercise the
    per-row `variant_id` -> joined `variant_ids` collapse behavior.
    """
    return dict(
        group_id=group_id,
        variant_id=variant_id if variant_id is not None else group_id,
        source_id='src',
        sample_id=sample_id,
        platform='pacbio',
        method=method,
        chromosome_1=chromosome_1 if chromosome_1 is not None else chromosome,
        position_1=position_1,
        strand_1='+',
        operation_1='D',
        chromosome_2=chromosome_2 if chromosome_2 is not None else chromosome,
        position_2=position_2,
        strand_2='+',
        operation_2='U',
        sequence=sequence,
        quality=quality,
        filter=filter_value,
        variant_size=variant_size,
        variant_type=variant_type,
        position_1_read_count_reference_allele=20,
        position_2_read_count_reference_allele=20,
        position_1_read_count_alternate_allele=p1_alt,
        position_2_read_count_alternate_allele=p2_alt,
        position_1_read_count_total=20 + p1_alt,
        position_2_read_count_total=20 + p2_alt,
        position_1_alternate_allele_fraction=0.3,
        position_2_alternate_allele_fraction=0.3,
        info='',
        format='',
        read_ids_alternate_allele='',
    )


def test_collapse_single_group_basic():
    """One group, one tumor row -> one collapsed row with informational cols appended."""
    df = pd.DataFrame([_row(1, 'tumor1')])
    out = collapse(
        df_variants=df,
        tumor_sample_ids=['tumor1'],
        normal_sample_ids=['normal1'],
        min_methods=1,
    )
    assert len(out) == 1
    assert out['variant_id'].iloc[0] == 1
    assert out['num_methods'].iloc[0] == 1
    assert out['num_rows_in_group'].iloc[0] == 1
    # Informational columns are appended at the end.
    assert list(out.columns)[-2:] == ['num_methods', 'num_rows_in_group']


def test_collapse_normal_contamination_drops_group():
    """Group with any normal-row alt-allele evidence is dropped entirely."""
    df = pd.DataFrame([
        # Group 1: clean (normal has 0 alt reads)
        _row(1, 'tumor1', p1_alt=8, p2_alt=8),
        _row(1, 'normal1', p1_alt=0, p2_alt=0),
        # Group 2: contaminated (normal has alt reads at position_1)
        _row(2, 'tumor1', p1_alt=8, p2_alt=8),
        _row(2, 'normal1', p1_alt=2, p2_alt=0),
        # Group 3: contaminated (normal has alt reads at position_2 only)
        _row(3, 'tumor1', p1_alt=8, p2_alt=8),
        _row(3, 'normal1', p1_alt=0, p2_alt=3),
    ])
    out = collapse(
        df_variants=df,
        tumor_sample_ids=['tumor1'],
        normal_sample_ids=['normal1'],
        min_methods=1,
    )
    assert sorted(out['variant_id'].tolist()) == [1]


def test_collapse_min_methods_filter():
    """Groups whose unique-method count is below the threshold are dropped."""
    df = pd.DataFrame([
        # Group 1: 1 method only
        _row(1, 'tumor1', method='savana'),
        # Group 2: 2 methods
        _row(2, 'tumor1', method='savana'),
        _row(2, 'tumor1', method='manta'),
        # Group 3: 3 methods
        _row(3, 'tumor1', method='savana'),
        _row(3, 'tumor1', method='manta'),
        _row(3, 'tumor1', method='delly'),
    ])
    out_min1 = collapse(df, ['tumor1'], ['normal1'], min_methods=1)
    out_min2 = collapse(df, ['tumor1'], ['normal1'], min_methods=2)
    out_min3 = collapse(df, ['tumor1'], ['normal1'], min_methods=3)
    assert sorted(out_min1['variant_id'].tolist()) == [1, 2, 3]
    assert sorted(out_min2['variant_id'].tolist()) == [2, 3]
    assert sorted(out_min3['variant_id'].tolist()) == [3]


def test_collapse_position_median_from_tumor_rows():
    """position_1 and position_2 should be the median across tumor rows, rounded to int.

    Importantly, normal-row positions must NOT influence the median.
    """
    df = pd.DataFrame([
        _row(1, 'tumor1', method='savana', position_1=1000, position_2=2000),
        _row(1, 'tumor1', method='manta',  position_1=1010, position_2=2020),
        _row(1, 'tumor1', method='delly',  position_1=1020, position_2=2040),
        # Normal row at very different positions — should be ignored for median.
        _row(1, 'normal1', method='savana', position_1=9999, position_2=9999, p1_alt=0, p2_alt=0),
    ])
    out = collapse(df, ['tumor1'], ['normal1'], min_methods=1)
    assert len(out) == 1
    assert out['position_1'].iloc[0] == 1010   # median(1000, 1010, 1020)
    assert out['position_2'].iloc[0] == 2020   # median(2000, 2020, 2040)


def test_collapse_numeric_max_from_tumor_rows():
    """Numeric stats (quality, read counts) collapse to max across tumor rows."""
    df = pd.DataFrame([
        _row(1, 'tumor1', method='savana', quality=10.0, p1_alt=5, p2_alt=4),
        _row(1, 'tumor1', method='manta',  quality=80.0, p1_alt=12, p2_alt=15),
        _row(1, 'tumor1', method='delly',  quality=50.0, p1_alt=8, p2_alt=10),
    ])
    out = collapse(df, ['tumor1'], ['normal1'], min_methods=1)
    assert out['quality'].iloc[0] == 80.0
    assert out['position_1_read_count_alternate_allele'].iloc[0] == 12
    assert out['position_2_read_count_alternate_allele'].iloc[0] == 15


def test_collapse_joined_columns_preserve_all_rows():
    """Joined columns dedup with first-occurrence order and include normal rows.

    This is what gives the user "all original information preserved" — the
    joined sample_id and method columns let you recover which samples and
    callers contributed to the group.
    """
    df = pd.DataFrame([
        _row(1, 'tumor1', method='savana'),
        _row(1, 'tumor1', method='manta'),
        _row(1, 'normal1', method='savana', p1_alt=0, p2_alt=0),
        _row(1, 'normal1', method='manta', p1_alt=0, p2_alt=0),
    ])
    out = collapse(df, ['tumor1'], ['normal1'], min_methods=1)
    # Both samples preserved, deduped.
    assert set(out['sample_id'].iloc[0].split('|')) == {'tumor1', 'normal1'}
    # Both methods preserved, deduped.
    assert set(out['method'].iloc[0].split('|')) == {'savana', 'manta'}
    # Counts reflect ALL rows in the group (tumor + normal), not just tumor.
    assert out['num_methods'].iloc[0] == 2
    assert out['num_rows_in_group'].iloc[0] == 4


def test_collapse_group_with_no_tumor_rows_is_dropped():
    """If a group has no rows matching --tumor-sample-ids, it cannot be represented and is dropped."""
    df = pd.DataFrame([
        # Group 1: tumor + normal -> kept
        _row(1, 'tumor1'),
        _row(1, 'normal1', p1_alt=0, p2_alt=0),
        # Group 2: only normal (no tumor row at all) -> dropped
        _row(2, 'normal1', p1_alt=0, p2_alt=0),
    ])
    out = collapse(df, ['tumor1'], ['normal1'], min_methods=1)
    assert out['variant_id'].tolist() == [1]


def test_collapse_group_with_no_normal_rows_is_kept():
    """No normal rows in a group means no contamination evidence — group survives."""
    df = pd.DataFrame([
        _row(1, 'tumor1'),
    ])
    out = collapse(df, ['tumor1'], ['normal1'], min_methods=1)
    assert len(out) == 1
    assert out['variant_id'].iloc[0] == 1


def test_collapse_no_normal_sample_ids_skips_contamination_filter():
    """Tumor-only callsets: omitting normal_sample_ids must not drop groups.

    A merged TSV may have no normal sample at all (e.g. tumor-only cohort).
    In that case the user should be able to call collapse() with normal_sample_ids
    as None or an empty list, and the contamination filter is simply skipped.
    Crucially, even tumor rows with > 0 alt reads must NOT be misread as
    "contamination" — only rows whose sample_id is in normal_sample_ids count.
    """
    df = pd.DataFrame([
        _row(1, 'tumor1', p1_alt=8, p2_alt=8),
        _row(2, 'tumor1', p1_alt=12, p2_alt=15),
    ])
    # Default (None) normal_sample_ids.
    out_default = collapse(df_variants=df, tumor_sample_ids=['tumor1'], min_methods=1)
    assert sorted(out_default['variant_id'].tolist()) == [1, 2]
    # Explicit empty list.
    out_empty = collapse(df, ['tumor1'], [], min_methods=1)
    assert sorted(out_empty['variant_id'].tolist()) == [1, 2]
    # Explicit None.
    out_none = collapse(df, ['tumor1'], None, min_methods=1)
    assert sorted(out_none['variant_id'].tolist()) == [1, 2]


def test_collapse_tumor_and_normal_overlap_raises():
    """A sample ID cannot be both tumor and normal."""
    df = pd.DataFrame([_row(1, 'tumor1')])
    with pytest.raises(ValueError, match='both'):
        collapse(
            df_variants=df,
            tumor_sample_ids=['tumor1', 'shared'],
            normal_sample_ids=['normal1', 'shared'],
            min_methods=1,
        )


def test_collapse_missing_group_id_raises():
    """Input must already be a merged TSV with a `group_id` column."""
    df = pd.DataFrame([_row(1, 'tumor1')]).drop(columns=['group_id'])
    with pytest.raises(ValueError, match='group_id'):
        collapse(
            df_variants=df,
            tumor_sample_ids=['tumor1'],
            normal_sample_ids=['normal1'],
            min_methods=1,
        )


def test_collapse_empty_result_preserves_columns():
    """When everything is filtered out, the output has the expected schema, just empty."""
    df = pd.DataFrame([
        _row(1, 'tumor1', method='savana'),
    ])
    # min_methods=2 will drop the only group (it has 1 method)
    out = collapse(df, ['tumor1'], ['normal1'], min_methods=2)
    assert len(out) == 0
    # Schema is preserved: original columns (renamed) + appended informational
    # columns. The group identifier is exposed as `variant_id`; the original
    # per-row `variant_id` becomes `variant_ids`; no `group_id` survives.
    assert 'variant_id' in out.columns
    assert 'variant_ids' in out.columns
    assert 'group_id' not in out.columns
    assert 'num_methods' in out.columns
    assert 'num_rows_in_group' in out.columns


def test_collapse_sequence_kept_from_first_row_not_joined():
    """`sequence` collapses to the representative (first tumor) row's value only.

    Pipe-joining sequences from multiple callers produced a multi-sequence cell
    (e.g. "ACGT|TTGG") that is not itself a valid sequence and broke downstream
    loaders. Collapse must keep just the first row's sequence.
    """
    df = pd.DataFrame([
        _row(1, 'tumor1', method='savana', sequence='ACGT'),
        _row(1, 'tumor1', method='manta',  sequence='TTGG'),
        _row(1, 'tumor1', method='delly',  sequence='CCAA'),
    ])
    out = collapse(df, ['tumor1'], min_methods=1)
    assert len(out) == 1
    seq = out['sequence'].iloc[0]
    assert seq == 'ACGT'          # first row only
    assert '|' not in str(seq)    # never concatenated


def test_collapse_filter_kept_from_first_row_not_joined():
    """`filter` collapses to the representative (first tumor) row's value only.

    Pipe-joining filter statuses across callers produced cells like "PASS|."
    that are not a meaningful single filter value and broke downstream loaders.
    Collapse must keep just the first row's filter.
    """
    df = pd.DataFrame([
        _row(1, 'tumor1', method='savana', filter_value='PASS'),
        _row(1, 'tumor1', method='manta',  filter_value='.'),
        _row(1, 'tumor1', method='delly',  filter_value='LowQual'),
    ])
    out = collapse(df, ['tumor1'], min_methods=1)
    assert len(out) == 1
    filt = out['filter'].iloc[0]
    assert filt == 'PASS'         # first row only
    assert '|' not in str(filt)   # never concatenated


def test_collapse_output_renames_group_id_and_joins_variant_ids():
    """Collapse output exposes the group identifier as `variant_id` and
    pipe-joins the original per-row `variant_id` values into `variant_ids`.

    Input is a `vstol merge` output: each row carries the group's `group_id`
    plus its own original `variant_id`. After collapse there is one row per
    group whose `variant_id` is the former group identifier and whose
    `variant_ids` lists the originating variant ids (deduped, first-occurrence
    order). The old `group_id` column must not remain.
    """
    df = pd.DataFrame([
        _row(1, 'tumor1', method='savana', variant_id=101),
        _row(1, 'tumor1', method='manta',  variant_id=202),
        # Duplicate original id within the group collapses to a single token.
        _row(1, 'tumor1', method='delly',  variant_id=101),
    ])
    out = collapse(df, ['tumor1'], min_methods=1)
    assert len(out) == 1
    r = out.iloc[0]
    # The group identifier is now the variant_id.
    assert r['variant_id'] == 1
    # Original per-row variant_ids are joined (deduped, in first-occurrence order).
    assert r['variant_ids'] == '101|202'
    # The pre-rename column name is gone.
    assert 'group_id' not in out.columns


def test_collapse_translocation_mixed_orientation_keeps_positions_paired():
    """Regression: a translocation grouped from BOTH breakend orders must keep
    each chromosome paired with its own position.

    A group can hold the same event as chr16<->chr5 and chr5<->chr16 (e.g.
    lumpy mate records, or two callers ordering breakends differently).
    Medianing position_1 across all rows averaged a chr16 coordinate with a
    chr5 one, leaving chromosome_1 (from the representative row) paired with a
    position on neither chromosome — "a position outside chr16's size". Each
    position's median must be taken only over rows on its representative
    chromosome.
    """
    df = pd.DataFrame([
        _row(1, 'tumor1', method='lumpy', variant_type='TRA', variant_size=-1,
             chromosome_1='chr16', position_1=1000,
             chromosome_2='chr5', position_2=180_000_000),
        _row(1, 'tumor1', method='manta', variant_type='TRA', variant_size=-1,
             chromosome_1='chr5', position_1=180_000_000,
             chromosome_2='chr16', position_2=1000),
    ])
    out = collapse(df, ['tumor1'], min_methods=1)
    assert len(out) == 1
    r = out.iloc[0]
    # Representative orientation is the first row: chr16 <-> chr5.
    assert r['chromosome_1'] == 'chr16' and r['position_1'] == 1000
    assert r['chromosome_2'] == 'chr5' and r['position_2'] == 180_000_000


def test_collapse_translocation_consistent_orientation_medians_jitter():
    """With a single consistent orientation, the per-chromosome median still
    smooths caller-specific breakpoint jitter."""
    df = pd.DataFrame([
        _row(1, 'tumor1', method='lumpy', variant_type='TRA', variant_size=-1,
             chromosome_1='chr16', position_1=1000,
             chromosome_2='chr5', position_2=5_000_000),
        _row(1, 'tumor1', method='manta', variant_type='TRA', variant_size=-1,
             chromosome_1='chr16', position_1=1002,
             chromosome_2='chr5', position_2=5_000_004),
        _row(1, 'tumor1', method='delly', variant_type='TRA', variant_size=-1,
             chromosome_1='chr16', position_1=1004,
             chromosome_2='chr5', position_2=5_000_006),
    ])
    out = collapse(df, ['tumor1'], min_methods=1)
    assert len(out) == 1
    r = out.iloc[0]
    assert r['chromosome_1'] == 'chr16' and r['position_1'] == 1002   # median 1000/1002/1004
    assert r['chromosome_2'] == 'chr5' and r['position_2'] == 5_000_004  # median 5M/5M+4/5M+6
