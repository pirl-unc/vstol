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
The purpose of this python3 script is to implement common functions related
to handling VCF files.
"""


import gzip
import pandas as pd
from typing import Dict
from ..constants import VariantCallingMethods


def read_vcf_file(
        vcf_file: str,
        low_memory=True,
        memory_map=False
) -> pd.DataFrame:
    """
    Read a VCF file and return a Pandas DataFrame.

    Parameters:
        vcf_file    :   VCF file.
        low_memory  :   Low memory (default: True).
        memory_map  :   Map memory (default: False).

    Returns:
        Pandas DataFrame
    """
    vcf_names = []
    is_gzipped = False
    if vcf_file.endswith(".gz"):
        is_gzipped = True
        with gzip.open(vcf_file, 'rt') as f:
            for line in f:
                if line.startswith("#CHROM"):
                    vcf_names = line.split('\t')
                    break
    else:
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith("#CHROM"):
                    vcf_names = line.split('\t')
                    break

    vcf_names = [i.replace('\n', '') for i in vcf_names]
    vcf_names = ['CHROM' if i == '#CHROM' else i for i in vcf_names]
    if is_gzipped:
        return pd.read_csv(vcf_file,
                           compression='gzip',
                           comment='#',
                           sep='\s+',
                           header=None,
                           low_memory=low_memory,
                           memory_map=memory_map,
                           names=vcf_names)
    else:
        return pd.read_csv(vcf_file,
                           comment='#',
                           sep='\s+',
                           header=None,
                           low_memory=low_memory,
                           memory_map=memory_map,
                           names=vcf_names)


def get_attribute_types(variant_calling_method: str) -> Dict:
    """
    Get the attribute type dictionary for a variant calling method.
    """
    if variant_calling_method == VariantCallingMethods.CUTESV:
        return VariantCallingMethods.AttributeTypes.CUTESV
    elif variant_calling_method == VariantCallingMethods.DEEPVARIANT:
        return VariantCallingMethods.AttributeTypes.DEEPVARIANT
    elif variant_calling_method == VariantCallingMethods.DELLY2_SOMATIC:
        return VariantCallingMethods.AttributeTypes.DELLY2_SOMATIC
    elif variant_calling_method == VariantCallingMethods.GATK4_MUTECT2:
        return VariantCallingMethods.AttributeTypes.GATK4_MUTECT2
    elif variant_calling_method == VariantCallingMethods.LUMPY_SOMATIC:
        return VariantCallingMethods.AttributeTypes.LUMPY_SOMATIC
    elif variant_calling_method == VariantCallingMethods.MANTA_SOMATIC:
        return VariantCallingMethods.AttributeTypes.MANTA_SOMATIC
    elif variant_calling_method == VariantCallingMethods.PBSV:
        return VariantCallingMethods.AttributeTypes.PBSV
    elif variant_calling_method == VariantCallingMethods.SAVANA:
        return VariantCallingMethods.AttributeTypes.SAVANA
    elif variant_calling_method == VariantCallingMethods.SEVERUS:
        return VariantCallingMethods.AttributeTypes.SEVERUS
    elif variant_calling_method == VariantCallingMethods.SNIFFLES2:
        return VariantCallingMethods.AttributeTypes.SNIFFLES2
    elif variant_calling_method == VariantCallingMethods.STRELKA2_SOMATIC:
        return VariantCallingMethods.AttributeTypes.STRELKA2_SOMATIC
    elif variant_calling_method == VariantCallingMethods.SVIM:
        return VariantCallingMethods.AttributeTypes.SVIM
    elif variant_calling_method == VariantCallingMethods.SVISIONPRO:
        return VariantCallingMethods.AttributeTypes.SVISIONPRO
    else:
        raise Exception('Unsupported variant calling method: %s')


def write_vcf_file(
        variants_list: 'VariantsList',
        output_vcf_file: str,
        gzip: bool = False
):
    """
    Write a VariantsList object as a VCF file.

    Parameters:
        variants_list       :   VariantsList object.
        output_vcf_file     :   Output VCF file.
        gzip                :   If True, gzip the VCF file.
    """
    # Step 1. Check whether the VariantsList is ready to be written as a VCF file
    # Make sure there is only one VariantCall object per Variant and that there
    # is only one sample ID.
    sample_ids = set()
    for variant in variants_list.variants:
        if variant.num_variant_calls != 1:
            raise Exception('Please make sure all variants have exactly one '
                            'variant call. Variant %s has %i variant calls.'
                            % (variant.id, variant.num_variant_calls))
        for variant_call in variant.variant_calls:
            sample_ids.add(variant_call.sample_id)
    if len(sample_ids) != 1:
        raise Exception('Please make sure there is only one sample ID. '
                        'Supplied variants list has %i sample IDs.'
                        % (len(sample_ids)))
    sample_id = list(sample_ids)[0]

    # Step 2. Write VCF
    if gzip:
        if output_vcf_file.endswith(".gz") == False:
            output_vcf_file = output_vcf_file + '.gz'
        f = gzip.open(output_vcf_file, 'wt')
    else:
        f = open(output_vcf_file, 'w')
    f.write('##fileformat=VCFv4.2\n')
    f.write('##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Variant with precise breakpoints">\n')
    f.write('##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Variant with imprecise breakpoints">\n')
    f.write('##INFO=<ID=SIZE,Number=1,Type=Integer,Description="Variant size">\n')
    f.write('##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type">\n')
    f.write('##INFO=<ID=CHR2,Number=1,Type=String,Description="Mate chromsome for BND SVs">\n')
    f.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variant">\n')
    f.write('##INFO=<ID=RNAMES,Number=.,Type=String,Description="Names of supporting reads">\n')
    f.write('##INFO=<ID=METHOD,Number=.,Type=String,Description="Variant calling method">\n')
    f.write('##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">\n')
    f.write('##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">\n')
    f.write('##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">\n')
    f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % sample_id)
    for variant in variants_list.variants:
        variant_call = variant.variant_calls[0]
        curr_chrom = variant_call.chromosome_1
        curr_pos = str(variant_call.position_1)
        curr_id = variant_call.id
        curr_ref = variant_call.reference_allele
        curr_alt = variant_call.alternate_allele
        curr_qual = '.' if variant_call.quality_score == -1 else str(variant_call.quality_score)
        curr_filter = variant_call.filter
        info = []
        if variant_call.precise:
            info.append('PRECISE')
        else:
            info.append('IMPREICSE')
        if variant_call.variant_size != -1:
            info.append('SIZE=%i' % variant_call.variant_size)
        info.append('TYPE=%s' % variant_call.variant_type)
        info.append('CHR2=%s' % variant_call.chromosome_2)
        info.append('END=%i' % variant_call.position_2)
        if len(variant_call.alternate_allele_read_ids) > 0:
            info.append('RNAMES=%s' % ','.join(list(variant_call.alternate_allele_read_ids)))
        info.append('METHOD=%s' % variant_call.variant_calling_method)
        curr_info = ';'.join(info)
        format = 'DR:DV:VAF'
        dr = '.' if variant_call.reference_allele_read_count == -1 else variant_call.reference_allele_read_count
        dv = '.' if variant_call.alternate_allele_read_count == -1 else variant_call.alternate_allele_read_count
        vaf = '.' if variant_call.alternate_allele_fraction == -1.0 else variant_call.alternate_allele_fraction
        curr_value = '%s:%s:%s' % (str(dr), str(dv), str(vaf))
        f.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'
                % (curr_chrom,
                   curr_pos,
                   curr_id,
                   curr_ref,
                   curr_alt,
                   curr_qual,
                   curr_filter,
                   curr_info,
                   format,
                   curr_value))
    f.close()
