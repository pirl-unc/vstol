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
        with gzip.open(vcf_file, 'r') as f:
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
                           delim_whitespace=True,
                           header=None,
                           low_memory=low_memory,
                           memory_map=memory_map,
                           names=vcf_names)
    else:
        return pd.read_csv(vcf_file,
                           comment='#',
                           delim_whitespace=True,
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
    elif variant_calling_method == VariantCallingMethods.PBSV:
        return VariantCallingMethods.AttributeTypes.PBSV
    elif variant_calling_method == VariantCallingMethods.SNIFFLES2:
        return VariantCallingMethods.AttributeTypes.SNIFFLES2
    elif variant_calling_method == VariantCallingMethods.STRELKA2_SOMATIC:
        return VariantCallingMethods.AttributeTypes.STRELKA2_SOMATIC
    elif variant_calling_method == VariantCallingMethods.SVIM:
        return VariantCallingMethods.AttributeTypes.SVIM
    else:
        raise Exception('Unsupported variant calling method: %s')


def write_vcf_file(
        vcf_file: str,
        variants_list: 'VariantsList'
):
    """
    Write a VariantsList object as a VCF file.

    Parameters:
        vcf_file        :   VCF file.
        variants_list   :   VariantsList.
    """
    # Step 1. Write header
    with open(vcf_file, 'w') as f:
        f.write('##fileformat=VCFv4.2\n')
        f.write('##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">\n')
        f.write('##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">\n')
        f.write('##FORMAT=<ID=VC,Number=R,Type=String,Description="Variant calling methods">\n')

