from vstolib.main import visualize
from .data import get_data_path


def test_visualize_cutesv():
    tsv_file = get_data_path(name='hg002_cutesv.tsv')
    pdf_file = get_data_path(name='hg002_cutesv.pdf')
    visualize(
        rscript_path='Rscript',
        tsv_file=tsv_file,
        output_pdf_file=pdf_file
    )

def test_visualize_deepvariant():
    tsv_file = get_data_path(name='hg002_deepvariant.tsv')
    pdf_file = get_data_path(name='hg002_deepvariant.pdf')
    visualize(
        rscript_path='Rscript',
        tsv_file=tsv_file,
        output_pdf_file=pdf_file
    )

def test_visualize_delly2():
    tsv_file = get_data_path(name='hg002_hg001_delly2.tsv')
    pdf_file = get_data_path(name='hg002_hg001_delly2.pdf')
    visualize(
        rscript_path='Rscript',
        tsv_file=tsv_file,
        output_pdf_file=pdf_file
    )

def test_visualize_manta():
    tsv_file = get_data_path(name='sample001_manta_somatic.tsv')
    pdf_file = get_data_path(name='sample001_manta_somatic.pdf')
    visualize(
        rscript_path='Rscript',
        tsv_file=tsv_file,
        output_pdf_file=pdf_file
    )

def test_visualize_mutect2():
    tsv_file = get_data_path(name='hg002_hg001_gatk4_mutect2.tsv')
    pdf_file = get_data_path(name='hg002_hg001_gatk4_mutect2.pdf')
    visualize(
        rscript_path='Rscript',
        tsv_file=tsv_file,
        output_pdf_file=pdf_file
    )

def test_visualize_lumpy():
    tsv_file = get_data_path(name='hg002_hg001_lumpy.tsv')
    pdf_file = get_data_path(name='hg002_hg001_lumpy.pdf')
    visualize(
        rscript_path='Rscript',
        tsv_file=tsv_file,
        output_pdf_file=pdf_file
    )

def test_visualize_pbsv():
    tsv_file = get_data_path(name='hg002_pbsv.tsv')
    pdf_file = get_data_path(name='hg002_pbsv.pdf')
    visualize(
        rscript_path='Rscript',
        tsv_file=tsv_file,
        output_pdf_file=pdf_file
    )

def test_visualize_savana():
    tsv_file = get_data_path(name='sample001_savana_sv_breakpoints.tsv')
    pdf_file = get_data_path(name='sample001_savana_sv_breakpoints.pdf')
    visualize(
        rscript_path='Rscript',
        tsv_file=tsv_file,
        output_pdf_file=pdf_file
    )

def test_visualize_severus():
    tsv_file = get_data_path(name='sample001_severus_somatic.tsv')
    pdf_file = get_data_path(name='sample001_severus_somatic.pdf')
    visualize(
        rscript_path='Rscript',
        tsv_file=tsv_file,
        output_pdf_file=pdf_file
    )

def test_visualize_sniffles2():
    tsv_file = get_data_path(name='hg002_sniffles2.tsv')
    pdf_file = get_data_path(name='hg002_sniffles2.pdf')
    visualize(
        rscript_path='Rscript',
        tsv_file=tsv_file,
        output_pdf_file=pdf_file
    )

def test_visualize_strelka2():
    tsv_file = get_data_path(name='hg002_hg001_strelka2_indels.tsv')
    pdf_file = get_data_path(name='hg002_hg001_strelka2_indels.pdf')
    visualize(
        rscript_path='Rscript',
        tsv_file=tsv_file,
        output_pdf_file=pdf_file
    )
    tsv_file = get_data_path(name='hg002_hg001_strelka2_snvs.tsv')
    pdf_file = get_data_path(name='hg002_hg001_strelka2_snvs.pdf')
    visualize(
        rscript_path='Rscript',
        tsv_file=tsv_file,
        output_pdf_file=pdf_file
    )

def test_visualize_svim():
    tsv_file = get_data_path(name='hg002_svim.tsv')
    pdf_file = get_data_path(name='hg002_svim.pdf')
    visualize(
        rscript_path='Rscript',
        tsv_file=tsv_file,
        output_pdf_file=pdf_file
    )

def test_visualize_svisionpro():
    tsv_file = get_data_path(name='sample001.svision_pro_v1.8.s3.tsv')
    pdf_file = get_data_path(name='sample001.svision_pro_v1.8.s3.pdf')
    visualize(
        rscript_path='Rscript',
        tsv_file=tsv_file,
        output_pdf_file=pdf_file
    )