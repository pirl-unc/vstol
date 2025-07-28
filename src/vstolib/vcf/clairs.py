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
The purpose of this python3 script is to implement a parser for ClairS VCF files.
"""


import pandas as pd
from .common import get_template_dict
from ..constants import *
from ..logging import get_logger
from ..utilities import retrieve_from_dict


logger = get_logger(__name__)


def parse_clairs_vcf(
        df_vcf: pd.DataFrame,
        source_id: str,
        platform: str
) -> pd.DataFrame:
    """
    Parse a ClairS VCF DataFrame and return a Pandas DataFrame.

    Parameters:
        df_vcf      :   Pandas DataFrame of rows from a ClairS VCF file.

    Returns:
        Pandas DataFrame with the following columns:
        'id'
        'source_id'
        'sample_id'
        'platform'
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
        'read_count_reference_allele'
        'read_count_alternate_allele'
        'read_count_total'
        'alternate_allele_fraction'
        'info'
        'format'
    """
    # Step 1. Get a template dictionary
    data = get_template_dict()

    # Step 2. Get all sample IDs
    sample_ids = df_vcf.columns.values.tolist()[9:]

    # Step 3. Record variant calls
    variant_id = 1
    for row in df_vcf.to_dict('records'):
        # Get VCF column values
        chromosome = retrieve_from_dict(dct=row, key='CHROM', default_value='', type=str)
        position = retrieve_from_dict(dct=row, key='POS', default_value='', type=int)
        reference_allele = retrieve_from_dict(dct=row, key='REF', default_value='', type=str)
        alternate_allele = retrieve_from_dict(dct=row, key='ALT', default_value='', type=str)
        quality = retrieve_from_dict(dct=row, key='QUAL', default_value=-1.0, type=float)
        filter = retrieve_from_dict(dct=row, key='FILTER', default_value='', type=str)
        info = str(row['INFO']).split(';')
        format = str(row['FORMAT']).split(':')

        # Get INFO
        info_dict = {}
        for curr_info in info:
            if '=' in curr_info:
                curr_info_elements = curr_info.split('=')
                curr_key = curr_info_elements[0]
                info_dict[curr_key] = str(curr_info_elements[1])
            else:
                info_dict[curr_info] = True

        for sample_id in sample_ids:
            # Get FORMAT
            format_dict = {}
            format_curr_sample = str(row[sample_id]).split(':')
            for curr_format in format:
                curr_key = curr_format
                format_dict[curr_key] = str(format_curr_sample[format.index(curr_format)])

            # Store meta data
            data['id'].append(variant_id)
            data['source_id'].append(source_id)
            data['sample_id'].append(sample_id)
            data['platform'].append(platform)
            data['method'].append(VariantCallingMethod.CLAIRS)
            data['quality'].append(quality)
            data['filter'].append(filter)

            # Store variant information
            data['chromosome_1'].append(chromosome)
            data['chromosome_2'].append(chromosome)
            data['strand_1'].append(Strand.UNKNOWN)
            data['strand_2'].append(Strand.UNKNOWN)
            data['operation_1'].append(OperationType.DOWNSTREAM)
            data['operation_2'].append(OperationType.UPSTREAM)
            if len(reference_allele) == len(alternate_allele):
                # SNV or MNV
                data['position_1'].append(position - 1)
                data['position_2'].append(position + 1)
                data['sequence'].append(alternate_allele)
                data['size'].append(len(alternate_allele))
                if len(alternate_allele) == 1:
                    data['variant_type'].append(VariantType.SINGLE_NUCLEOTIDE_VARIANT)
                else:
                    data['variant_type'].append(VariantType.MULTI_NUCLEOTIDE_VARIANT)
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

            # Store FORMAT
            format_key_value_pairs = []
            for key, val in format_dict.items():
                format_key_value_pairs.append('%s=%s' % (str(key), str(val)))
            data['format'].append(';'.join(format_key_value_pairs))

            # Store optional values
            if 'AD' in format_dict:
                position_1_read_count_reference_allele = int(format_dict['AD'].split(',')[0])
                position_2_read_count_reference_allele = int(format_dict['AD'].split(',')[0])
                position_1_read_count_alternate_allele = int(format_dict['AD'].split(',')[1])
                position_2_read_count_alternate_allele = int(format_dict['AD'].split(',')[1])
                position_1_read_count_total = position_1_read_count_reference_allele + position_1_read_count_alternate_allele
                position_2_read_count_total = position_2_read_count_reference_allele + position_2_read_count_alternate_allele
            else:
                position_1_read_count_reference_allele = -1
                position_2_read_count_reference_allele = -1
                position_1_read_count_alternate_allele = -1
                position_2_read_count_alternate_allele = -1
                position_1_read_count_total = -1
                position_2_read_count_total = -1
            if 'AF' in format_dict:
                position_1_alternate_allele_fraction = float(format_dict['AF'])
                position_2_alternate_allele_fraction = float(format_dict['AF'])
            else:
                position_1_alternate_allele_fraction = -1.0
                position_2_alternate_allele_fraction = -1.0
            data['position_1_read_count_reference_allele'].append(position_1_read_count_reference_allele)
            data['position_2_read_count_reference_allele'].append(position_2_read_count_reference_allele)
            data['position_1_read_count_alternate_allele'].append(position_1_read_count_alternate_allele)
            data['position_2_read_count_alternate_allele'].append(position_2_read_count_alternate_allele)
            data['position_1_read_count_total'].append(position_1_read_count_total)
            data['position_2_read_count_total'].append(position_2_read_count_total)
            data['position_1_alternate_allele_fraction'].append(position_1_alternate_allele_fraction)
            data['position_2_alternate_allele_fraction'].append(position_2_alternate_allele_fraction)
            data['read_ids_alternate_allele'].append('')
        variant_id += 1

    df_variants = pd.DataFrame(data)
    logger.info('Converted %i ClairS variants.' % len(df_variants['id'].unique()))
    return df_variants
