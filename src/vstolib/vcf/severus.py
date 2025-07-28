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
The purpose of this python3 script is to implement a parser for Severus VCF files.
"""


from .common import *
from ..constants import *
from ..logging import get_logger
from ..utilities import retrieve_from_dict


logger = get_logger(__name__)


def parse_severus_vcf(
        df_vcf: pd.DataFrame,
        source_id: str,
        platform: str
) -> pd.DataFrame:
    """
    Parse a Severus VCF DataFrame and return a Pandas DataFrame.

    Parameters:
        df_vcf      :   Pandas DataFrame of rows from a Severus VCF file.

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
    included_ids = set()
    for row in df_vcf.to_dict('records'):
        # Get VCF column values
        chromosome = retrieve_from_dict(dct=row, key='CHROM', default_value='', type=str)
        position = retrieve_from_dict(dct=row, key='POS', default_value='', type=int)
        reference_allele = retrieve_from_dict(dct=row, key='REF', default_value='', type=str)
        id = retrieve_from_dict(dct=row, key='ID', default_value='', type=str)
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

        if info_dict['SVTYPE'] == 'BND':
            if info_dict['MATE_ID'] in included_ids:
                continue

        included_ids.add(id)

        for i, sample_id in enumerate(sample_ids):
            # Get FORMAT
            format_dict = {}
            format_curr_sample = str(row[sample_id]).split(':')
            for curr_format in format:
                curr_key = curr_format
                format_dict[curr_key] = str(format_curr_sample[format.index(curr_format)])

            # Read counts and allele fractions
            dr = int(format_dict['DR']) if format_dict.get('DR', '.') != '.' else -1
            dv = int(format_dict['DV']) if format_dict.get('DV', '.') != '.' else -1
            af = float(format_dict.get('VAF', '.').split(',')[i]) if format_dict.get('VAF', '.').split(',')[i] != '.' else -1.0

            def append_variant(
                    chromosome_1,
                    position_1,
                    strand_1,
                    operation_1,
                    chromosome_2,
                    position_2,
                    strand_2,
                    operation_2,
                    sequence,
                    size,
                    variant_type
            ):
                data['id'].append(variant_id)
                data['source_id'].append(source_id)
                data['sample_id'].append(sample_id)
                data['platform'].append(platform)
                data['method'].append(VariantCallingMethod.SEVERUS)
                data['quality'].append(quality)
                data['filter'].append(filter)
                data['chromosome_1'].append(chromosome_1)
                data['position_1'].append(position_1)
                data['strand_1'].append(strand_1)
                data['operation_1'].append(operation_1)
                data['chromosome_2'].append(chromosome_2)
                data['position_2'].append(position_2)
                data['strand_2'].append(strand_2)
                data['operation_2'].append(operation_2)
                data['sequence'].append(sequence)
                data['size'].append(size)
                data['variant_type'].append(variant_type)
                data['info'].append(';'.join('%s=%s' % (k, v) for k, v in info_dict.items()))
                data['format'].append(';'.join('%s=%s' % (k, v) for k, v in format_dict.items()))
                data['position_1_read_count_reference_allele'].append(dr)
                data['position_2_read_count_reference_allele'].append(dr)
                data['position_1_read_count_alternate_allele'].append(dv)
                data['position_2_read_count_alternate_allele'].append(dv)
                data['position_1_read_count_total'].append(dr + dv if dr != -1 and dv != -1 else -1)
                data['position_2_read_count_total'].append(dr + dv if dr != -1 and dv != -1 else -1)
                data['position_1_alternate_allele_fraction'].append(af)
                data['position_2_alternate_allele_fraction'].append(af)
                data['read_ids_alternate_allele'].append('')

            sv_type = info_dict['SVTYPE']
            if sv_type == 'DEL':
                deletion_start = position + 1
                deletion_end = int(info_dict['END'])
                append_variant(
                    chromosome_1=chromosome,
                    position_1=position - 1,
                    strand_1=Strand.FORWARD,
                    operation_1=OperationType.DOWNSTREAM,
                    chromosome_2=chromosome,
                    position_2=deletion_end + 1,
                    strand_2=Strand.FORWARD,
                    operation_2=OperationType.UPSTREAM,
                    sequence='',
                    size=deletion_end - deletion_start + 1,
                    variant_type=VariantType.DELETION
                )
            elif sv_type == 'INS':
                inserted_sequence = alternate_allele[1:]
                append_variant(
                    chromosome_1=chromosome,
                    position_1=position,
                    strand_1=Strand.FORWARD,
                    operation_1=OperationType.DOWNSTREAM,
                    chromosome_2=chromosome,
                    position_2=position + 1,
                    strand_2=Strand.FORWARD,
                    operation_2=OperationType.UPSTREAM,
                    sequence=inserted_sequence,
                    size=len(inserted_sequence),
                    variant_type=VariantType.INSERTION
                )
            elif sv_type in ['BND', 'TRA']:
                alt_breakpoints = parse_alt_breakpoints(
                    alternate_allele=alternate_allele,
                    chromosome_1=chromosome,
                    position_1=position
                )
                variant_type = VariantType.BREAKPOINT if alt_breakpoints['chromosome_1'] == alt_breakpoints['chromosome_2'] else VariantType.TRANSLOCATION
                size = abs(alt_breakpoints['position_1'] - alt_breakpoints['position_2']) if variant_type == VariantType.BREAKPOINT else -1
                append_variant(
                    chromosome_1=alt_breakpoints['chromosome_1'],
                    position_1=alt_breakpoints['position_1'],
                    strand_1=alt_breakpoints['strand_1'],
                    operation_1=alt_breakpoints['operation_1'],
                    chromosome_2=alt_breakpoints['chromosome_2'],
                    position_2=alt_breakpoints['position_2'],
                    strand_2=alt_breakpoints['strand_2'],
                    operation_2=alt_breakpoints['operation_2'],
                    sequence=alt_breakpoints['sequence'],
                    size=size,
                    variant_type=variant_type
                )
            elif sv_type == 'DUP':
                duplication_end = int(info_dict['END'])
                size = int(info_dict['SVLEN'])
                append_variant(
                    chromosome_1=chromosome,
                    position_1=position,
                    strand_1=Strand.UNKNOWN,
                    operation_1=OperationType.UPSTREAM,
                    chromosome_2=chromosome,
                    position_2=duplication_end,
                    strand_2=Strand.UNKNOWN,
                    operation_2=OperationType.DOWNSTREAM,
                    sequence='',
                    size=size,
                    variant_type=VariantType.DUPLICATION
                )
            elif sv_type == 'INV':
                if info_dict['STRANDS'] == '++':
                    inversion_start = position + 1
                    inversion_end = inversion_start + int(info_dict['SVLEN']) - 1
                    size = abs(inversion_end - inversion_start) + 1
                else:
                    inversion_start = position
                    inversion_end = inversion_start + int(info_dict['SVLEN']) - 1
                    size = abs(inversion_end - inversion_start) + 1

                # Head-to-head inversion breakpoint
                append_variant(
                    chromosome_1=chromosome,
                    position_1=inversion_start - 1,
                    strand_1=Strand.FORWARD,
                    operation_1=OperationType.DOWNSTREAM,
                    chromosome_2=chromosome,
                    position_2=inversion_end,
                    strand_2=Strand.REVERSE,
                    operation_2=OperationType.DOWNSTREAM,
                    sequence='',
                    size=size,
                    variant_type=VariantType.INVERSION
                )

                # Tail-to-tail inversion breakpoint
                append_variant(
                    chromosome_1=chromosome,
                    position_1=inversion_start,
                    strand_1=Strand.REVERSE,
                    operation_1=OperationType.UPSTREAM,
                    chromosome_2=chromosome,
                    position_2=inversion_end + 1,
                    strand_2=Strand.FORWARD,
                    operation_2=OperationType.UPSTREAM,
                    sequence='',
                    size=size,
                    variant_type=VariantType.INVERSION
                )
            else:
                raise Exception('Unexpected SV type: %s' % sv_type)

        variant_id += 1

    df_variants = pd.DataFrame(data)
    logger.info('Converted %i Severus variants.' % len(df_variants['id'].unique()))
    return df_variants

