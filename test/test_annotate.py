from vstolib.ensembl import Ensembl
from vstolib.gencode import Gencode
from vstolib.main import annotate
from .data import get_data_path


def test_annotate_sniffles2_ensembl(sniffles2_variants_list):
    ensembl = Ensembl(release=95, species='human')
    variants_list_annotated = annotate(
        variants_list=sniffles2_variants_list,
        annotator=ensembl
    )
    print(variants_list_annotated.size)

def test_annotate_sniffles2_gencode(sniffles2_variants_list):
    gtf_file = get_data_path(name='gencode.v41.annotations.gtf')
    gencode = Gencode(
        gtf_file=gtf_file,
        version='v41',
        species='human',
    )
    variants_list_annotated = annotate(
        variants_list=sniffles2_variants_list,
        annotator=gencode
    )
    print(variants_list_annotated.size)
