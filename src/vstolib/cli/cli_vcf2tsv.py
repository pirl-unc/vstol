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
The purpose of this python3 script is to create parser
and run 'vcf2tsv' command.
"""


import argparse
from ..constants import *
from ..logging import get_logger
from ..main import vcf2tsv
from ..vcf.common import read_vcf_file


logger = get_logger(__name__)


def add_cli_vcf2tsv_arg_parser(
        sub_parsers: argparse._SubParsersAction
) -> argparse._SubParsersAction:
    """
    Add 'vcf2tsv' parser.
    """
    parser = sub_parsers.add_parser('vcf2tsv', help='Convert a VCF file to a TSV file.')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--vcf-file", '-i',
        dest="vcf_file",
        type=str,
        required=True,
        help="Input VCF file."
    )
    parser_required.add_argument(
        "--method", '-m',
        dest="method",
        type=str,
        required=True,
        choices=VariantCallingMethod.ALL,
        help="Variant calling method. "
             "Allowed options: %s."
             % (', '.join(VariantCallingMethod.ALL))
    )
    parser_required.add_argument(
        "--platform", '-p',
        dest="platform",
        type=str,
        required=True,
        help="Sequencing platform."
    )
    parser_required.add_argument(
        "--source-id", '-s',
        dest="source_id",
        type=str,
        required=True,
        help="Source ID (e.g. patient ID or cell line ID)."
    )
    parser_required.add_argument(
        "--output-tsv-file", '-o',
        dest="output_tsv_file",
        type=str,
        required=True,
        help="Output TSV file."
    )

    parser.set_defaults(which='vcf2tsv')
    return sub_parsers


def run_cli_vcf2tsv_from_parsed_args(args: argparse.Namespace):
    """
    Run 'vcf2tsv' command using parameters from parsed arguments.

    Parameters:
        args    :   argparse.Namespace object with the following variables:
                    vcf_file
                    method
                    platform
                    source_id
                    output_tsv_file
    """
    df_vcf = read_vcf_file(vcf_file=args.vcf_file)
    df_variants = vcf2tsv(
        df_vcf=df_vcf,
        source_id=args.source_id,
        method=args.method,
        platform=args.platform
    )
    df_variants.to_csv(args.output_tsv_file, sep='\t', index=False)
