import pandas as pd
from vstolib.gencode import Gencode
from vstolib.main import annotate
from .data import get_data_path


def test_annotate_1():
    gencode = Gencode(
        gtf_file=get_data_path(name="gtf/gencode.v41.annotation.chr17-18.gtf.gz"),
        version='v41',
        species='human',
        levels=[1,2],
        types=['protein_coding']
    )

    data = {
        'id': [1],
        'chromosome_1': ['chr17'],
        'position_1': [7673800],
        'strand_1': ['*'],
        'chromosome_2': ['chr17'],
        'position_2': [7673801],
        'strand_2': ['*'],
    }

    df_variants = pd.DataFrame(data)

    df_annotated = annotate(
        df_variants=df_variants,
        annotator=gencode,
        num_processes=1
    )

    assert df_annotated['position_1_exon_id'].values[0] != ''
    assert df_annotated['position_2_exon_id'].values[0] != ''


def test_annotate_2():
    gencode = Gencode(
        gtf_file=get_data_path(name="gtf/gencode.v41.annotation.chr17-18.gtf.gz"),
        version='v41',
        species='human',
        levels=[1,2],
        types=['protein_coding']
    )

    data = {
        'id': [1],
        'chromosome_1': ['chr17'],
        'position_1': [7673800],
        'strand_1': ['-'],
        'chromosome_2': ['chr17'],
        'position_2': [7673801],
        'strand_2': ['-'],
    }

    df_variants = pd.DataFrame(data)

    df_annotated = annotate(
        df_variants=df_variants,
        annotator=gencode,
        num_processes=1
    )

    assert df_annotated['position_1_exon_id'].values[0] != ''
    assert df_annotated['position_2_exon_id'].values[0] != ''


def test_annotate_3():
    gencode = Gencode(
        gtf_file=get_data_path(name="gtf/gencode.v41.annotation.chr17-18.gtf.gz"),
        version='v41',
        species='human',
        levels=[1,2],
        types=['protein_coding']
    )

    data = {
        'id': [1],
        'chromosome_1': ['chr17'],
        'position_1': [7673800],
        'strand_1': ['+'],
        'chromosome_2': ['chr17'],
        'position_2': [7673801],
        'strand_2': ['+'],
    }

    df_variants = pd.DataFrame(data)

    df_annotated = annotate(
        df_variants=df_variants,
        annotator=gencode,
        num_processes=1
    )

    assert df_annotated['position_1_exon_id'].values[0] == ''
    assert df_annotated['position_2_exon_id'].values[0] == ''
