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
and run 'tsv2vcf' command.
"""


import argparse
from ..constants import *
from ..default import *
from ..logging import get_logger
from ..main import collapse
from ..utilities import str2bool
from ..variants_list import VariantsList
from ..vcf.common import write_vcf_file


logger = get_logger(__name__)


def add_cli_tsv2vcf_arg_parser(
        sub_parsers: argparse._SubParsersAction
) -> argparse._SubParsersAction:
    """
    Add 'tsv2vcf' parser.
    """
    parser = sub_parsers.add_parser('tsv2vcf', help='Convert a TSV file to a VCF file.')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--tsv-file", '-i',
        dest="tsv_file",
        type=str,
        required=True,
        help="Input TSV file."
    )
    parser_required.add_argument(
        "--sample-id", '-a',
        dest="sample_id",
        type=str,
        required=True,
        help="Sample ID to retain."
    )
    parser_required.add_argument(
        "--output-vcf-file", '-o',
        dest="output_vcf_file",
        type=str,
        required=True,
        help="Output VCF file."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--strategy",
        dest="strategy",
        type=str,
        required=False,
        default=CollapseStrategies.MAX_ALTERNATE_ALLELE_READ_COUNT,
        choices=CollapseStrategies.ALL,
        help="Collapse (summarization) strategy. "
             "Allowed options: %s. Default: %s"
             % (', '.join(CollapseStrategies.ALL),
                CollapseStrategies.MAX_ALTERNATE_ALLELE_READ_COUNT)
    )
    parser_optional.add_argument(
        "--gzip",
        dest="gzip",
        type=str2bool,
        required=False,
        default=GZIP,
        help="If 'yes', gzip the output TSV file (default: %s)."
             % GZIP
    )

    parser.set_defaults(which='tsv2vcf')
    return sub_parsers


def run_cli_tsv2vcf_from_parsed_args(args: argparse.Namespace):
    """
    Run 'tsv2vcf' command using parameters from parsed arguments.

    Parameters:
        args    :   argparse.Namespace object with the following variables:
                    tsv_file
                    sample_id
                    output_vcf_file
                    strategy
                    gzip
    """
    variants_list = VariantsList.read_tsv_file(tsv_file=args.tsv_file)
    variants_list_collapsed = collapse(
        variants_list=variants_list,
        sample_id=args.sample_id,
        strategy=args.strategy
    )
    write_vcf_file(
        variants_list=variants_list_collapsed,
        output_vcf_file=args.output_vcf_file,
        gzip=args.gzip
    )
