# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


"""
The purpose of this python3 script is to implement a parser for dbSNP VCF files.
"""


import pandas as pd
from .common import get_template_dict
from ..constants import OperationType, Strand, VariantCallingMethod, VariantType
from ..logging import get_logger
from ..utilities import retrieve_from_dict


logger = get_logger(__name__)


def parse_dbsnp_vcf(
        df_vcf: pd.DataFrame
) -> pd.DataFrame:
    """
    Parse a dBSNP VCF DataFrame and return a Pandas DataFrame.

    Parameters:
        df_vcf      :   Pandas DataFrame of rows from a ClairS VCF file.

    Returns:
        Pandas DataFrame with the following columns:
        'id'
        'method'
        'chromosome_1'
        'position_1'
        'strand_1'
        'operation_1'
        'chromosome_2'
        'position_2'
        'strand_2'
        'operation_2'
        'sequence'
        'quality'
        'filter'
        'size'
        'info'
    """
    # Step 1. Get a template dictionary
    data = get_template_dict()
    del data['source_id']
    del data['sample_id']
    del data['platform']
    del data['format']
    del data['read_count_reference_allele']
    del data['read_count_alternate_allele']
    del data['read_count_total']
    del data['alternate_allele_fraction']
    del data['read_ids_alternate_allele']

    # Step 2. Record variant calls
    for row in df_vcf.to_dict('records'):
        # Get VCF column values
        chromosome = retrieve_from_dict(dct=row, key='CHROM', default_value='', type=str)
        position = retrieve_from_dict(dct=row, key='POS', default_value='', type=int)
        variant_id = retrieve_from_dict(dct=row, key='ID', default_value='', type=str)
        reference_allele = retrieve_from_dict(dct=row, key='REF', default_value='', type=str)
        alternate_allele = retrieve_from_dict(dct=row, key='ALT', default_value='', type=str)
        quality = retrieve_from_dict(dct=row, key='QUAL', default_value=-1.0, type=float)
        filter = retrieve_from_dict(dct=row, key='FILTER', default_value='', type=str)
        info = str(row['INFO']).split(';')

        # Get INFO
        info_dict = {}
        for curr_info in info:
            if '=' in curr_info:
                curr_info_elements = curr_info.split('=')
                curr_key = curr_info_elements[0]
                info_dict[curr_key] = str(curr_info_elements[1])
            else:
                info_dict[curr_info] = True

        # Store meta data
        data['id'].append(variant_id)
        data['method'].append(VariantCallingMethod.DBSNP)
        data['quality'].append(quality)
        data['filter'].append(filter)

        # Store variant information
        data['chromosome_1'].append(chromosome)
        data['chromosome_2'].append(chromosome)
        data['strand_1'].append(Strand.UNKNOWN)
        data['strand_2'].append(Strand.UNKNOWN)
        data['operation_1'].append(OperationType.DOWNSTREAM)
        data['operation_2'].append(OperationType.UPSTREAM)
        if len(reference_allele) == 1 and len(alternate_allele) == 1:
            # SNV
            data['position_1'].append(position - 1)
            data['position_2'].append(position + 1)
            data['sequence'].append(alternate_allele)
            data['size'].append(1)
            data['variant_type'].append(VariantType.SINGLE_NUCLEOTIDE_VARIANT)
        elif (len(reference_allele) > len(alternate_allele)) and len(alternate_allele) == 1:
            # Deletion
            deletion_start = position + 1
            deletion_end = position + len(reference_allele) - 1
            data['position_1'].append(deletion_start - 1)
            data['position_2'].append(deletion_end + 1)
            data['sequence'].append('')
            data['size'].append(deletion_end - deletion_start + 1)
            data['variant_type'].append(VariantType.DELETION)
        elif (len(reference_allele) < len(alternate_allele)) and len(reference_allele) == 1:
            # Insertion
            data['position_1'].append(position)
            data['position_2'].append(position + 1)
            data['sequence'].append(alternate_allele[1:])
            data['size'].append(len(alternate_allele[1:]))
            data['variant_type'].append(VariantType.INSERTION)
        else:
            raise Exception('Unexpected case. REF: %s. ALT: %s.' % (reference_allele, alternate_allele))

        # Store INFO
        info_key_value_pairs = []
        for key, val in info_dict.items():
            info_key_value_pairs.append('%s=%s' % (str(key), str(val)))
        data['info'].append(';'.join(info_key_value_pairs))

    df_variants = pd.DataFrame(data)
    logger.info('Converted %i dbSNP variants.' % len(df_variants['id'].unique()))
    return df_variants

