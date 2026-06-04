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
The purpose of this python3 script is to implement a parser for SVABA VCF files.
"""


from .common import *
from ..constants import *
from ..logging import get_logger
from ..utilities import retrieve_from_dict


logger = get_logger(__name__)


def parse_svaba_vcf(
        df_vcf: pd.DataFrame,
        source_id: str,
        platform: str
) -> pd.DataFrame:
    """
    Parse a SVABA VCF DataFrame and return a Pandas DataFrame.
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
                curr_info_elements = curr_info.split('=', 1)
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
                format_dict[curr_key] = str(format_curr_sample[format.index(curr_format)]) if format.index(curr_format) < len(format_curr_sample) else '.'

            # Read counts and allele fractions
            num_alt_reads = int(format_dict['AD']) if format_dict.get('AD','.').split(',') != '.' else -1
            num_total_reads = int(format_dict['DP']) if format_dict.get('DP','.').split(',') != '.' else -1
            if num_alt_reads > 0 and num_total_reads > 0:
                num_ref_reads = num_total_reads - num_alt_reads
                af = float(num_alt_reads) / float(num_total_reads)
            else:
                num_ref_reads = -1
                af = -1.0
            try:
                read_names = info_dict['READNAMES']
            except:
                read_names = ''

            if 'SVTYPE' in info_dict:
                # Structural variants
                if info_dict['SVTYPE'] == 'BND':
                    if info_dict['MATEID'] in included_ids:
                        continue

                included_ids.add(id)

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
                    data['variant_id'].append(variant_id)
                    data['source_id'].append(source_id)
                    data['sample_id'].append(sample_id)
                    data['platform'].append(platform)
                    data['method'].append(VariantCallingMethod.SVABA)
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
                    data['variant_size'].append(size)
                    data['variant_type'].append(variant_type)
                    data['info'].append(';'.join('%s=%s' % (k, v) for k, v in info_dict.items()))
                    data['format'].append(';'.join('%s=%s' % (k, v) for k, v in format_dict.items()))
                    data['position_1_read_count_reference_allele'].append(num_ref_reads)
                    data['position_2_read_count_reference_allele'].append(num_ref_reads)
                    data['position_1_read_count_alternate_allele'].append(num_alt_reads)
                    data['position_2_read_count_alternate_allele'].append(num_alt_reads)
                    data['position_1_read_count_total'].append(num_total_reads)
                    data['position_2_read_count_total'].append(num_total_reads)
                    data['position_1_alternate_allele_fraction'].append(af)
                    data['position_2_alternate_allele_fraction'].append(af)
                    data['read_ids_alternate_allele'].append(read_names)

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
                    insertion_sequence = info_dict['SVINSSEQ']
                    append_variant(
                        chromosome_1=chromosome,
                        position_1=position,
                        strand_1=Strand.FORWARD,
                        operation_1=OperationType.DOWNSTREAM,
                        chromosome_2=chromosome,
                        position_2=position + 1,
                        strand_2=Strand.FORWARD,
                        operation_2=OperationType.UPSTREAM,
                        sequence=insertion_sequence,
                        size=len(insertion_sequence),
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
                    size = abs(int(info_dict['END']) - position)
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
                else:
                    raise Exception('Unexpected SV type: %s' % sv_type)
            else:
                # Insertions and deletions
                data['variant_id'].append(variant_id)
                data['source_id'].append(source_id)
                data['sample_id'].append(sample_id)
                data['platform'].append(platform)
                data['method'].append(VariantCallingMethod.SVABA)
                data['quality'].append(quality)
                data['filter'].append(filter)
                data['chromosome_1'].append(chromosome)
                data['chromosome_2'].append(chromosome)
                data['strand_1'].append(Strand.FORWARD)
                data['strand_2'].append(Strand.FORWARD)
                data['operation_1'].append(OperationType.DOWNSTREAM)
                data['operation_2'].append(OperationType.UPSTREAM)
                if len(reference_allele) == len(alternate_allele):
                    # SNV or MNV
                    data['position_1'].append(position - 1)
                    data['position_2'].append(position + 1)
                    data['sequence'].append(alternate_allele)
                    data['variant_size'].append(len(alternate_allele))
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
                    data['variant_size'].append(deletion_end - deletion_start + 1)
                    data['variant_type'].append(VariantType.DELETION)
                elif (len(reference_allele) < len(alternate_allele)) and len(reference_allele) == 1:
                    # Insertion
                    data['position_1'].append(position)
                    data['position_2'].append(position + 1)
                    data['sequence'].append(alternate_allele[1:])
                    data['variant_size'].append(len(alternate_allele[1:]))
                    data['variant_type'].append(VariantType.INSERTION)
                else:
                    raise Exception('Unexpected case. REF: %s. ALT: %s.' % (reference_allele, alternate_allele))

                # Store optional values
                data['info'].append(';'.join('%s=%s' % (k, v) for k, v in info_dict.items()))
                data['format'].append(';'.join('%s=%s' % (k, v) for k, v in format_dict.items()))
                data['position_1_read_count_reference_allele'].append(num_ref_reads)
                data['position_2_read_count_reference_allele'].append(num_ref_reads)
                data['position_1_read_count_alternate_allele'].append(num_alt_reads)
                data['position_2_read_count_alternate_allele'].append(num_alt_reads)
                data['position_1_read_count_total'].append(num_total_reads)
                data['position_2_read_count_total'].append(num_total_reads)
                data['position_1_alternate_allele_fraction'].append(af)
                data['position_2_alternate_allele_fraction'].append(af)
                data['read_ids_alternate_allele'].append(read_names)

        variant_id += 1

    df_variants = pd.DataFrame(data)
    logger.info('Converted %i SVABA variants.' % len(df_variants['variant_id'].unique()))
    return df_variants

