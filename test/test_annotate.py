from vstolib.ensembl import Ensembl
from vstolib.gencode import Gencode
from vstolib.main import annotate
from python.vstolib.constants import GenomicRegionTypes
from .data import get_data_path


def test_annotate_gencode():
    gencode = Gencode(
        gtf_file=get_data_path(name="gencode.v41.annotations.gtf"),
        version='v41',
        species='human',
        levels=[1,2],
        types=['protein_coding']
    )

    # 5prime_utr
    variation_annotations = gencode.annotate_position(
        chromosome='chr1',
        position=1324620
    )
    for v in variation_annotations:
        assert(v.region == GenomicRegionTypes.FIVE_PRIME_UTR)

    # Exonic
    variation_annotations = gencode.annotate_position(
        chromosome='chr1',
        position=1313820
    )
    for v in variation_annotations:
        assert(v.region == GenomicRegionTypes.EXONIC)

   # Intronic
    variation_annotations = gencode.annotate_position(
        chromosome='chr1',
        position=1313327
    )
    for v in variation_annotations:
        assert(v.region == GenomicRegionTypes.INTRONIC)

   # Intergenic
    variation_annotations = gencode.annotate_position(
        chromosome='chr1',
        position=1329900
    )
    for v in variation_annotations:
        assert(v.region == GenomicRegionTypes.INTERGENIC)


def test_annotate_ensembl():
    ensembl = Ensembl(release=95, species='human')

    # UTR
    variant_annotations = ensembl.annotate_position_using_pyensembl(
        chromosome='chr8',
        position=94641250
    )
    for v in variant_annotations:
        assert(v.region == GenomicRegionTypes.UNTRANSLATED_REGION)

    # Exonic
    variant_annotations = ensembl.annotate_position_using_pyensembl(
        chromosome='chr8',
        position=94641400
    )
    for v in variant_annotations:
        assert(v.region == GenomicRegionTypes.EXONIC)

    # Intronic
    variant_annotations = ensembl.annotate_position_using_pyensembl(
        chromosome='chr8',
        position=94641600
    )
    for v in variant_annotations:
        assert(v.region == GenomicRegionTypes.INTRONIC)

    # Intergenic
    variant_annotations = ensembl.annotate_position_using_pyensembl(
        chromosome='chr8',
        position=94640800
    )
    for v in variant_annotations:
        assert(v.region == GenomicRegionTypes.INTERGENIC)

    # Exonic
    variant_annotations = ensembl.annotate_position_using_pyensembl(
        chromosome='chr1',
        position=12030
    )
    for v in variant_annotations:
        assert(v.region == GenomicRegionTypes.EXONIC)

    # Intronic
    variant_annotations = ensembl.annotate_position_using_pyensembl(
        chromosome='chr1',
        position=12300
    )
    for v in variant_annotations:
        assert(v.region == GenomicRegionTypes.INTRONIC)


def test_annotate_sniffles2_ensembl(sniffles2_variants_list):
    ensembl = Ensembl(release=95, species='human')
    variants_list_annotated = annotate(
        variants_list=sniffles2_variants_list,
        annotator=ensembl,
        num_threads=2
    )
    print(variants_list_annotated.size)


def test_annotate_sniffles2_gencode(sniffles2_variants_list):
    gtf_file = get_data_path(name='gencode.v41.annotations.gtf')
    gencode = Gencode(
        gtf_file=gtf_file,
        version='v41',
        species='human',
        levels=[1],
        types=['protein_coding']
    )
    variants_list_annotated = annotate(
        variants_list=sniffles2_variants_list,
        annotator=gencode,
        num_threads=2
    )
    print(variants_list_annotated.size)
