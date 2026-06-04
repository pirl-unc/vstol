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
import re
from ..constants import OperationType, Strand


def get_template_dict() -> dict:
    data = {
        'variant_id': [],
        'source_id': [],
        'sample_id': [],
        'platform': [],
        'method': [],
        'chromosome_1': [],
        'position_1': [],
        'strand_1': [],
        'operation_1': [],
        'chromosome_2': [],
        'position_2': [],
        'strand_2': [],
        'operation_2': [],
        'sequence': [],
        'quality': [],
        'filter': [],
        'variant_size': [],
        'variant_type': [],
        'position_1_read_count_reference_allele': [],
        'position_2_read_count_reference_allele': [],
        'position_1_read_count_alternate_allele': [],
        'position_2_read_count_alternate_allele': [],
        'position_1_read_count_total': [],
        'position_2_read_count_total': [],
        'position_1_alternate_allele_fraction': [],
        'position_2_alternate_allele_fraction': [],
        'info': [],
        'format': [],
        'read_ids_alternate_allele': []
    }
    return data


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
    num_header_lines = 0
    is_gzipped = vcf_file.endswith(".gz") or vcf_file.endswith(".bgz")
    open_func = gzip.open if is_gzipped else open
    with open_func(vcf_file, 'rt') as f:
        for line in f:
            num_header_lines += 1
            if line.startswith("#CHROM"):
                vcf_names = line.split('\t')
                break

    vcf_names = [i.replace('\n', '') for i in vcf_names]
    vcf_names = ['CHROM' if i == '#CHROM' else i for i in vcf_names]

    # Do NOT pass comment='#' to read_csv: pandas treats '#' anywhere in a line as a
    # comment, so a '#' inside a data row (INFO/FORMAT/read-name fields) would silently
    # truncate the row and drop its remaining columns. Skip exactly the header block
    # (meta lines + the #CHROM line) via skiprows instead.
    read_csv_kwargs = dict(
        sep='\t',
        header=None,
        skiprows=num_header_lines,
        low_memory=low_memory,
        memory_map=memory_map,
        names=vcf_names,
    )
    if is_gzipped:
        return pd.read_csv(vcf_file, compression='gzip', **read_csv_kwargs)
    else:
        return pd.read_csv(vcf_file, **read_csv_kwargs)


def parse_alt_breakpoints(
        alternate_allele: str,
        chromosome_1: str,
        position_1: int
) -> dict:
    """
    Get breakpoint orientation.

    Returns:
        dictionary = {
          'chromosome_1': <value>,
          'position_1': <value>,
          'strand_1': <value>,
          'operation_1': <value>,
          'chromosome_2': <value>,
          'position_2': <value>,
          'strand_2': <value>,
          'operation_2': <value>,
          'sequence': <value>
        }
    """
    # Single breakend (no brackets, dot at beginning or end)
    # Examples: "G." or ".A"
    has_brackets = ("[" in alternate_allele) or ("]" in alternate_allele)
    if not has_brackets and "." in alternate_allele:
        if alternate_allele.endswith("."):
            full = alternate_allele[:-1]
            sequence = full[1:] if len(full) > 1 else ""
            return {
                "chromosome_1": chromosome_1,
                "position_1": int(position_1),
                "strand_1": Strand.FORWARD,
                "operation_1": OperationType.DOWNSTREAM,
                "chromosome_2": "", # Unknown
                "position_2": -1,   # Unknown
                "strand_2": Strand.UNKNOWN,
                "operation_2": OperationType.UNKNOWN,
                "sequence": sequence
            }
        if alternate_allele.startswith("."):
            full = alternate_allele[1:]
            sequence = full[:-1] if len(full) > 1 else ""
            return {
                "chromosome_1": "",
                "position_1": -1,
                "strand_1": Strand.UNKNOWN,
                "operation_1": OperationType.UNKNOWN,
                "chromosome_2": chromosome_1,
                "position_2": int(position_1),
                "strand_2": Strand.FORWARD,
                "operation_2": OperationType.UPSTREAM,
                "sequence": sequence
            }
        raise Exception("Dot in an unexpected place: %s" % alternate_allele)

    # Standard BND with brackets
    pattern = re.compile(
        r"""
        (?P<prefix>[ACGTNacgtn]*)              # optional prefix bases
        (?P<left_bracket>[\[\]])               # [ or ]
        (?P<chromosome_2>[^:\[\]]+)            # anything up to the colon (no colon/brackets)
        :(?P<position_2>\d+)                   # position
        (?P<right_bracket>[\[\]])              # [ or ]
        (?P<suffix>[ACGTNacgtn]*)              # optional suffix bases
        """,
        re.VERBOSE
    )
    match = pattern.match(alternate_allele)
    groups = match.groupdict()

    if groups['prefix'] != '' and groups['left_bracket'] == '[' and groups['right_bracket'] == '[' and groups['suffix'] == '':
        # t[p[ piece extending to the right of p is joined after t
        # example: A[1:123456[
        chromosome_1 = chromosome_1
        position_1 = position_1
        strand_1 = Strand.FORWARD
        operation_1 = OperationType.DOWNSTREAM
        chromosome_2 = groups['chromosome_2']
        position_2 = groups['position_2']
        strand_2 = Strand.FORWARD
        operation_2 = OperationType.UPSTREAM
        if len(groups['prefix']) > 1:
            sequence = groups['prefix'][1:]
        else:
            sequence = ''
    elif groups['prefix'] != '' and groups['left_bracket'] == ']' and groups['right_bracket'] == ']' and groups['suffix'] == '':
        # t]p] reverse comp piece extending left of p is joined after t
        # example: A]1:123456]
        chromosome_1 = chromosome_1
        position_1 = position_1
        strand_1 = Strand.FORWARD
        operation_1 = OperationType.DOWNSTREAM
        chromosome_2 = groups['chromosome_2']
        position_2 = groups['position_2']
        strand_2 = Strand.REVERSE
        operation_2 = OperationType.DOWNSTREAM
        if len(groups['prefix']) > 1:
            sequence = groups['prefix'][1:]
        else:
            sequence = ''
    elif groups['prefix'] == '' and groups['left_bracket'] == ']' and groups['right_bracket'] == ']' and groups['suffix'] != '':
        # ]p]t piece extending to the left of p is joined before t
        # example: ]1:123456]A
        chromosome_2 = chromosome_1
        position_2 = position_1
        strand_2 = Strand.FORWARD
        operation_2 = OperationType.UPSTREAM
        chromosome_1 = groups['chromosome_2']
        position_1 = groups['position_2']
        strand_1 = Strand.FORWARD
        operation_1 = OperationType.DOWNSTREAM
        if len(groups['suffix']) > 1:
            sequence = groups['suffix'][:-1]
        else:
            sequence = ''
    elif groups['prefix'] == '' and groups['left_bracket'] == '[' and groups['right_bracket'] == '[' and groups['suffix'] != '':
        # [p[t  reverse comp piece extending right of p is joined before t
        # example: [1:123456[A
        chromosome_1 = chromosome_1
        position_1 = position_1
        strand_1 = Strand.FORWARD
        operation_1 = OperationType.UPSTREAM
        chromosome_2 = groups['chromosome_2']
        position_2 = groups['position_2']
        strand_2 = Strand.REVERSE
        operation_2 = OperationType.UPSTREAM
        if len(groups['suffix']) > 1:
            sequence = groups['suffix'][:-1]
        else:
            sequence = ''
    else:
        chromosome_1 = chromosome_1
        position_1 = position_1
        strand_1 = Strand.FORWARD
        operation_1 = OperationType.UNKNOWN
        chromosome_2 = groups['chromosome_2']
        position_2 = groups['position_2']
        strand_2 = Strand.REVERSE
        operation_2 = OperationType.UNKNOWN
        if len(groups['suffix']) > 1:
            sequence = groups['suffix'][:-1]
        else:
            sequence = ''

    return {
        'chromosome_1': chromosome_1,
        'position_1': int(position_1),
        'strand_1': strand_1,
        'operation_1': operation_1,
        'chromosome_2': chromosome_2,
        'position_2': int(position_2),
        'strand_2': strand_2,
        'operation_2': operation_2,
        'sequence': sequence
    }
