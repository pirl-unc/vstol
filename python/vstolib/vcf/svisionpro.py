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
The purpose of this python3 script is to implement a parser for SVision-pro VCF files.
"""


import pandas as pd
import re
from collections import OrderedDict
from typing import Dict
from ..constants import NucleicAcidTypes, VariantCallingMethods, VariantTypes
from ..logging import get_logger
from ..utilities import get_typed_value, retrieve_from_dict, retrieve_from_list
from ..variant import Variant
from ..variant_call import VariantCall
from ..variants_list import VariantsList


logger = get_logger(__name__)


def parse_svisionpro_callset(
        df_vcf: pd.DataFrame,
        sequencing_platform: str,
        source_id: str,
        case_id: str,
        control_id: str
) -> VariantsList:
    """
    Parses a SVision-pro DataFrame and returns a VariantsList object.

    Parameters:
        df_vcf                  :   DataFrame of rows from a SVision-pro VCF file.
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
        for idx, sample_id in enumerate(sample_ids):
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
                    curr_type = VariantCallingMethods.AttributeTypes.SVISIONPRO[curr_key]
                    attributes[curr_key] = get_typed_value(value=curr_info_elements[1], default_value='', type=curr_type)
                else:
                    if curr_info == 'PRECISE':
                        attributes['PRECISE'] = True
                    elif curr_info == 'IMPRECISE':
                        attributes['PRECISE'] = False
                    else:
                        attributes[curr_info] = True

            # Step 3. Extract FORMAT
            format = str(row['FORMAT']).split(':')
            curr_sample = str(row[sample_id]).split(':')
            for curr_format in format:
                curr_key = curr_format
                curr_type = VariantCallingMethods.AttributeTypes.SVISIONPRO[curr_key]
                attributes[curr_key] = retrieve_from_list(lst=curr_sample,
                                                          index=format.index(curr_format),
                                                          default_value='',
                                                          type=curr_type)

            # Step 4. Update variables
            # Update the following variables:
            # precise
            # variant_size
            # variant_type
            # chromosome_2
            # position_2
            # total_read_count
            # reference_allele_read_count
            # alternate_allele_read_count
            # alternate_allele_fraction
            # alternate_allele_read_ids
            if 'PRECISE' in attributes.keys():
                if attributes['PRECISE']:
                    precise = 'yes'
                else:
                    precise = 'no'
            if 'END' in attributes.keys():
                position_2 = attributes['END']
            if 'SVLEN' in attributes.keys():
                variant_size = abs(attributes['SVLEN'])
            if 'SVTYPE' in attributes.keys():
                if alternate_allele == 'CSV':
                    variant_type = 'CSV'
                else:
                    if attributes['SVTYPE'] == 'tDUP' or attributes['SVTYPE'] == 'dDUP':
                        variant_type = VariantTypes.DUPLICATION
                    else:
                        variant_type = attributes['SVTYPE']
            if 'DR' in attributes.keys():
                reference_allele_read_count = get_typed_value(value=attributes['DR'], default_value=-1, type=int)
            if 'DV' in attributes.keys():
                alternate_allele_read_count = get_typed_value(value=attributes['DV'], default_value=-1, type=int)
            if reference_allele_read_count != -1 and alternate_allele_read_count != -1:
                total_read_count = reference_allele_read_count + alternate_allele_read_count
            if alternate_allele_read_count != -1 and total_read_count != -1:
                alternate_allele_fraction = float(alternate_allele_read_count) / float(total_read_count)
            if 'RNAMES' in attributes.keys():
                alternate_allele_read_ids = attributes['RNAMES'].split(',')

            # Update variant_size for 'BND'
            if (variant_type == VariantTypes.BREAKPOINT or variant_type == VariantTypes.TRANSLOCATION) and \
                    (chromosome_1 == chromosome_2):
                variant_size = abs(position_1 - position_2) + 1

            # Append case variant call to variants list
            if idx == 0:
                sample_id = case_id
            else:
                sample_id = control_id
            variant_call_id = '%s_%s_%s_%i_%s_%s:%i_%s:%i' % (
                sample_id,
                NucleicAcidTypes.DNA,
                VariantCallingMethods.SVISIONPRO,
                curr_variant_call_idx,
                variant_type,
                chromosome_1,
                position_1,
                chromosome_2,
                position_2
            )
            curr_variant_call_idx += 1
            variant_call = VariantCall(
                id=variant_call_id,
                source_id=source_id,
                sample_id=sample_id,
                nucleic_acid=NucleicAcidTypes.DNA,
                variant_calling_method=VariantCallingMethods.SVISIONPRO,
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
        curr_variant_idx += 1

    variants_list = VariantsList(variants=list(variants.values()))
    logger.info('%i variants and %i variant calls in the returning VariantsList.' %
                (len(variants_list.variant_ids), len(variants_list.variant_call_ids)))
    return variants_list
