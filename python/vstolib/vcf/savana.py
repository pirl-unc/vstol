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
The purpose of this python3 script is to implement a parser for Savana VCF files.
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


def parse_savana_callset(
        df_vcf: pd.DataFrame,
        sequencing_platform: str,
        source_id: str,
        case_id: str,
        control_id: str
) -> VariantsList:
    """
    Parses a Savana DataFrame and returns a VariantsList object.

    Parameters:
        df_vcf                  :   DataFrame of rows from a Savana VCF file.
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
    included_mate_ids = set()
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
            precise = False
            tumor_total_read_count = -1
            tumor_reference_allele_read_count = -1
            tumor_alternate_allele_read_count = -1
            tumor_alternate_allele_fraction = -1.0
            normal_total_read_count = -1
            normal_reference_allele_read_count = -1
            normal_alternate_allele_read_count = -1
            normal_alternate_allele_fraction = -1.0
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
                    curr_type = VariantCallingMethods.AttributeTypes.SAVANA[curr_key]
                    attributes[curr_key] = get_typed_value(value=curr_info_elements[1], default_value='', type=curr_type)
                else:
                    attributes[curr_info] = True

            # Step 3. Extract FORMAT
            format = str(row['FORMAT']).split(':')
            curr_sample = str(row[sample_id]).split(':')
            for curr_format in format:
                curr_key = curr_format
                curr_type = VariantCallingMethods.AttributeTypes.SAVANA[curr_key]
                attributes[curr_key] = retrieve_from_list(lst=curr_sample,
                                                          index=format.index(curr_format),
                                                          default_value='',
                                                          type=curr_type)

            # Step 4. Update variables
            # Update the following variables:
            # variant_size
            # variant_type
            # chromosome_2
            # position_2
            # tumor_total_read_count
            # tumor_reference_allele_read_count
            # tumor_alternate_allele_read_count
            # tumor_alternate_allele_fraction
            # normal_total_read_count
            # normal_reference_allele_read_count
            # normal_alternate_allele_read_count
            # normal_alternate_allele_fraction
            if 'SVLEN' in attributes.keys():
                variant_size = abs(int(attributes['SVLEN']))
            if 'SVTYPE' in attributes.keys():
                variant_type = attributes['SVTYPE']
            if 'TUMOUR_DP' in attributes.keys():
                if ',' in attributes['TUMOUR_DP']: # SVTYPE=BND
                    tumor_total_read_count = int(float(attributes['TUMOUR_DP'].split(',')[0]))
                else :
                    tumor_total_read_count = int(float(attributes['TUMOUR_DP']))
            if 'NORMAL_DP' in attributes.keys():
                if ',' in attributes['NORMAL_DP']: # SVTYPE=BND
                    normal_total_read_count = int(float(attributes['NORMAL_DP'].split(',')[0]))
                else :
                    normal_total_read_count = int(float(attributes['NORMAL_DP']))
            if 'TUMOUR_SUPPORT' in attributes.keys():
                tumor_alternate_allele_read_count = int(attributes['TUMOUR_SUPPORT'])
            if 'NORMAL_SUPPORT' in attributes.keys():
                normal_alternate_allele_read_count = int(attributes['NORMAL_SUPPORT'])
            if tumor_total_read_count != -1 and tumor_alternate_allele_read_count != -1:
                tumor_reference_allele_read_count = tumor_total_read_count - tumor_alternate_allele_read_count
            if normal_total_read_count != -1 and normal_alternate_allele_read_count != -1:
                normal_reference_allele_read_count = normal_total_read_count - normal_alternate_allele_read_count
            if tumor_alternate_allele_read_count != -1 and tumor_total_read_count > 0:
                tumor_alternate_allele_fraction = float(tumor_alternate_allele_read_count) / float(tumor_total_read_count)
            if normal_alternate_allele_read_count != -1 and normal_total_read_count > 0:
                normal_alternate_allele_fraction = float(normal_alternate_allele_read_count) / float(normal_total_read_count)

            # Update chromosome_2 and position_2 for 'BND'
            if variant_type == VariantTypes.BREAKPOINT or variant_type == VariantTypes.TRANSLOCATION:
                pattern = re.compile(r'(chr\S+):(\d+)')
                matches = pattern.findall(str(row['ALT']))
                chromosome_2 = str(matches[0][0])
                position_2 = int(matches[0][1])
            else:
                position_2 = position_1

            # Update variant_size for 'BND'
            if (variant_type == VariantTypes.BREAKPOINT or variant_type == VariantTypes.TRANSLOCATION) and \
                    (chromosome_1 == chromosome_2):
                variant_size = abs(position_1 - position_2) + 1

            # Check if variant_call ID has been included
            if variant_type == VariantTypes.BREAKPOINT or variant_type == VariantTypes.TRANSLOCATION:
                if attributes['ID'] in included_mate_ids:
                    continue
                included_mate_ids.add(attributes['MATEID'])

            # Append case variant call to variants list
            case_variant_call_id = '%s_%s_%s_%i_%s_%s:%i_%s:%i' % (
                case_id,
                NucleicAcidTypes.DNA,
                VariantCallingMethods.SAVANA,
                curr_variant_call_idx,
                variant_type,
                chromosome_1,
                position_1,
                chromosome_2,
                position_2
            )
            curr_variant_call_idx += 1
            control_variant_call_id = '%s_%s_%s_%i_%s_%s:%i_%s:%i' % (
                control_id,
                NucleicAcidTypes.DNA,
                VariantCallingMethods.SAVANA,
                curr_variant_call_idx,
                variant_type,
                chromosome_1,
                position_1,
                chromosome_2,
                position_2
            )
            case_variant_call = VariantCall(
                id=case_variant_call_id,
                source_id=source_id,
                sample_id=case_id,
                nucleic_acid=NucleicAcidTypes.DNA,
                variant_calling_method=VariantCallingMethods.SAVANA,
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
                total_read_count=tumor_total_read_count,
                reference_allele_read_count=tumor_reference_allele_read_count,
                alternate_allele_read_count=tumor_alternate_allele_read_count,
                alternate_allele_fraction=tumor_alternate_allele_fraction,
                alternate_allele_read_ids=set(alternate_allele_read_ids),
                attributes=attributes
            )
            control_variant_call = VariantCall(
                id=control_variant_call_id,
                source_id=source_id,
                sample_id=control_id,
                nucleic_acid=NucleicAcidTypes.DNA,
                variant_calling_method=VariantCallingMethods.SAVANA,
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
                total_read_count=normal_total_read_count,
                reference_allele_read_count=normal_reference_allele_read_count,
                alternate_allele_read_count=normal_alternate_allele_read_count,
                alternate_allele_fraction=normal_alternate_allele_fraction,
                alternate_allele_read_ids=set(alternate_allele_read_ids),
                attributes=attributes
            )
            variant_id = str(curr_variant_idx)
            if variant_id not in variants:
                variants[variant_id] = Variant(id=variant_id)
            variants[variant_id].add_variant_call(variant_call=case_variant_call)
            variants[variant_id].add_variant_call(variant_call=control_variant_call)
        curr_variant_idx += 1

    variants_list = VariantsList(variants=list(variants.values()))
    logger.info('%i variants and %i variant calls in the returning VariantsList.' %
                (len(variants_list.variant_ids), len(variants_list.variant_call_ids)))
    return variants_list
