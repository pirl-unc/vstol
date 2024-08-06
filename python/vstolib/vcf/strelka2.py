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
The purpose of this python3 script is to implement a parser for Strelka2 VCF files.
"""


import pandas as pd
from collections import OrderedDict
from typing import Dict
from ..constants import NucleicAcidTypes, VariantCallingMethods, VariantTypes
from ..logging import get_logger
from ..utilities import get_typed_value, retrieve_from_dict, retrieve_from_list
from ..variant import Variant
from ..variant_call import VariantCall
from ..variants_list import VariantsList


logger = get_logger(__name__)


def parse_strelka2_somatic_callset(
        df_vcf: pd.DataFrame,
        sequencing_platform: str,
        source_id: str,
        case_id: str,
        control_id: str
) -> VariantsList:
    """
    Parses a Strelka2 DataFrame and returns VariantsList object.

    Parameters:
        df_vcf                  :   DataFrame of rows from a Strelka2 (germline mode) VCF file.
        sequencing_platform     :   Sequencing platform.
        source_id               :   Source ID.
        case_id                 :   Case ID.
        control_id              :   Control ID.

    Returns:
        VariantsList
    """
    variants: Dict[str, Variant] = {}
    sample_ids = df_vcf.columns.values.tolist()[9:]
    curr_variant_call_idx = 1
    curr_variant_idx = 1
    for row in df_vcf.to_dict('records'):
        for sample_id in sample_ids:
            # Step 1. Initialize values
            chromosome_1 = retrieve_from_dict(dct=row, key='CHROM', default_value='', type=str)
            chromosome_2 = retrieve_from_dict(dct=row, key='CHROM', default_value='', type=str)
            position_1 = retrieve_from_dict(dct=row, key='POS', default_value=-1, type=int)
            position_2 = -1
            reference_allele = retrieve_from_dict(dct=row, key='REF', default_value='', type=str)
            alternate_allele = retrieve_from_dict(dct=row, key='ALT', default_value='', type=str)
            filter = retrieve_from_dict(dct=row, key='FILTER', default_value='', type=str)
            quality_score = retrieve_from_dict(dct=row, key='QUAL', default_value=-1.0, type=float)
            precise = ''
            total_read_count = -1
            reference_allele_read_count = -1
            alternate_allele_read_count = -1
            alternate_allele_fraction = -1.0
            variant_type = ''
            variant_subtype = ''
            variant_sequences = []
            variant_size = -1
            alternate_allele_read_ids = []
            attributes = OrderedDict()
            attributes['ID'] = retrieve_from_dict(dct=row, key='ID', default_value='', type=str)

            # Step 2. Extract INFO
            info = str(row['INFO']).split(';')
            for curr_info in info:
                if '=' in curr_info:
                    curr_info_elements = curr_info.split('=')
                    curr_key = curr_info_elements[0]
                    curr_type = VariantCallingMethods.AttributeTypes.STRELKA2_SOMATIC[curr_key]
                    attributes[curr_key] = get_typed_value(value=curr_info_elements[1],
                                                           default_value='',
                                                           type=curr_type)
                else:
                    if curr_info == 'PRECISE':
                        attributes['PRECISE'] = True
                    elif curr_info == 'IMPRECISE':
                        attributes['PRECISE'] = False
                    attributes[curr_info] = True

            # Step 3. Extract FORMAT
            format = str(row['FORMAT']).split(':')
            curr_sample = str(row[sample_id]).split(':')
            for curr_format in format:
                curr_key = curr_format
                curr_type = VariantCallingMethods.AttributeTypes.STRELKA2_SOMATIC[curr_key]
                attributes[curr_key] = retrieve_from_list(lst=curr_sample,
                                                          index=format.index(curr_format),
                                                          default_value='',
                                                          type=curr_type)

            # Step 4. Update variables
            # Update the following variables:
            # precise
            # variant_type
            # variant_size
            # variant_sequences
            # position_1
            # position_2
            if 'PRECISE' in attributes.keys():
                if attributes['PRECISE']:
                    precise = 'yes'
                else:
                    precise = 'no'
            if len(reference_allele) == 1 and len(alternate_allele) == 1:
                variant_type = VariantTypes.SINGLE_NUCLEOTIDE_VARIANT
                variant_sequences.append(str(alternate_allele))
                variant_size = 1
                position_2 = position_1
            elif len(reference_allele) == 1 and len(alternate_allele) > 1:
                variant_type = VariantTypes.INSERTION
                variant_sequences.append(str(alternate_allele[1:]))
                variant_size = len(alternate_allele[1:])
                position_2 = position_1
            elif len(reference_allele) > 1 and len(alternate_allele) == 1:
                variant_type = VariantTypes.DELETION
                variant_sequences.append(str(reference_allele[1:]))
                variant_size = len(reference_allele[1:])
                position_1 = position_1 + 1
                position_2 = position_1 + variant_size - 1
            elif len(reference_allele) > 1 and len(alternate_allele) > 1:
                variant_type = VariantTypes.MULTI_NUCLEOTIDE_VARIANT
                variant_sequences.append(str(alternate_allele))
                variant_size = len(alternate_allele)
                position_2 = position_1 + len(variant_sequences[0]) - 1
            else:
                raise Exception(
                    'Unknown variant type. REF: %s. ALT: %s' %
                    (reference_allele, alternate_allele)
                )

            # Update the following variables:
            # total_read_count
            # reference_allele_read_count
            # alternate_allele_read_count
            if 'DP' in attributes.keys():
                total_read_count = attributes['DP']
            if 'AD' in attributes.keys():
                reference_allele_read_count = get_typed_value(value=attributes['AD'].split(',')[0],
                                                              default_value=-1,
                                                              type=int)
                alternate_allele_read_count = get_typed_value(value=attributes['AD'].split(',')[1],
                                                              default_value=-1,
                                                              type=int)

            if variant_type == VariantTypes.SINGLE_NUCLEOTIDE_VARIANT:
                if alternate_allele == 'A' or alternate_allele == 'a':
                    attr_key = 'AU'
                elif alternate_allele == 'C' or alternate_allele == 'c':
                    attr_key = 'CU'
                elif alternate_allele == 'G' or alternate_allele == 'g':
                    attr_key = 'GU'
                elif alternate_allele == 'T' or alternate_allele == 't':
                    attr_key = 'TU'
                else:
                    attr_key = ''
                if attr_key in attributes.keys():
                    alternate_allele_read_count = get_typed_value(value=attributes[attr_key].split(',')[0],
                                                                  default_value=-1,
                                                                  type=int)
            if variant_type == VariantTypes.INSERTION or VariantTypes.DELETION:
                if 'TIR' in attributes.keys():
                    alternate_allele_read_count = get_typed_value(value=attributes['TIR'].split(',')[0],
                                                                  default_value=-1,
                                                                  type=int)

            # Update the following variable if they are currently unknown but can be inferred
            # total_read_count
            # alternate_allele_fraction
            if reference_allele_read_count >= 0 and \
                alternate_allele_read_count >= 0 and \
                total_read_count == -1:
                total_read_count = alternate_allele_read_count + reference_allele_read_count
            if alternate_allele_read_count >= 0 and \
                total_read_count >= 0 and \
                reference_allele_read_count == -1:
                reference_allele_read_count = total_read_count - alternate_allele_read_count
            if alternate_allele_read_count >= 0 and \
                total_read_count > 0 and \
                alternate_allele_fraction == -1.0:
                alternate_allele_fraction = float(alternate_allele_read_count) / float(total_read_count)

            # Replace sample_id
            if sample_id == 'TUMOR':
                sample_id = case_id
            if sample_id == 'NORMAL':
                sample_id = control_id

            # Append variant call to variants list
            variant_call_id = '%s_%s_%s_%i_%s_%s:%i_%s:%i' % (
                sample_id,
                NucleicAcidTypes.DNA,
                VariantCallingMethods.STRELKA2_SOMATIC,
                curr_variant_call_idx,
                variant_type,
                chromosome_1,
                position_1,
                chromosome_2,
                position_2
            )
            variant_call = VariantCall(
                id=variant_call_id,
                source_id=source_id,
                sample_id=sample_id,
                nucleic_acid=NucleicAcidTypes.DNA,
                variant_calling_method=VariantCallingMethods.STRELKA2_SOMATIC,
                sequencing_platform=sequencing_platform,
                chromosome_1=chromosome_1,
                position_1=position_1,
                chromosome_2=chromosome_2,
                position_2=position_2,
                precise=precise,
                reference_allele=reference_allele,
                alternate_allele=alternate_allele,
                filter=filter,
                quality_score=quality_score,
                variant_type=variant_type,
                variant_subtype=variant_subtype,
                variant_size=variant_size,
                variant_sequences=set(variant_sequences),
                total_read_count=total_read_count,
                reference_allele_read_count=reference_allele_read_count,
                alternate_allele_read_count=alternate_allele_read_count,
                alternate_allele_fraction=alternate_allele_fraction,
                alternate_allele_read_ids=set(alternate_allele_read_ids),
                attributes=attributes
            )
            variant_id = str(curr_variant_idx)
            if variant_id not in variants:
                variants[variant_id] = Variant(id=variant_id)
            variants[variant_id].add_variant_call(variant_call=variant_call)
            curr_variant_call_idx += 1
        curr_variant_idx += 1

    variants_list = VariantsList(variants=list(variants.values()))
    logger.info('%i variants and %i variant calls in the returning VariantsList.' %
                (len(variants_list.variant_ids), len(variants_list.variant_call_ids)))
    return variants_list
