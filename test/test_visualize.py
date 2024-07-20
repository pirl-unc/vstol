from vstolib.main import visualize
from .data import get_data_path


def test_visualize(severus_variants_list):
    tsv_file = get_data_path(name='hg002_sniffles2.tsv')
    pdf_file = get_data_path(name='hg002_sniffles2.pdf')
    visualize(
        rscript_path='Rscript',
        tsv_file=tsv_file,
        output_pdf_file=pdf_file
    )