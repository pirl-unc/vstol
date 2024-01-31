from .conftest import *
from vstolib.constants import CollapseStrategies
from vstolib.main import collapse
from vstolib.variants_list import VariantsList


def test_collapse():
    tsv_file = get_data_path(name='hg002_merged_variants.tsv')
    variants_list = VariantsList.read_tsv_file(tsv_file=tsv_file)
    variants_list_collapsed = collapse(
        variants_list=variants_list,
        sample_id='HG002',
        strategy=CollapseStrategies.MAX_ALTERNATE_ALLELE_READ_COUNT
    )
    print(variants_list_collapsed.size)