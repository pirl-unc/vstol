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
and run 'visualize' command.
"""


import argparse
from ..constants import *
from ..logging import get_logger
from ..main import visualize
from ..vcf.common import read_vcf_file


logger = get_logger(__name__)


def add_cli_visuzlize_arg_parser(
        sub_parsers: argparse._SubParsersAction
) -> argparse._SubParsersAction:
    """
    Add 'visualize' parser.
    """
    parser = sub_parsers.add_parser('visualize', help='Visualizes a VSTOL TSV file.')
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
        "--output-pdf-file", '-o',
        dest="output_pdf_file",
        type=str,
        required=True,
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--rscript-path",
        dest="rscript_path",
        type=str,
        default='Rscript',
        required=False,
        help="Rscript path (default: Rscript)."
    )

    parser.set_defaults(which='visualize')
    return sub_parsers


def run_cli_visualize_from_parsed_args(args: argparse.Namespace):
    """
    Run 'visualize' command using parameters from parsed arguments.

    Parameters:
        args    :   argparse.Namespace object with the following variables:
                    tsv_file
                    output_pdf_file
    """
    visualize(
        rscript_path=args.rscript_path,
        tsv_file=args.tsv_file,
        output_pdf_file=args.output_pdf_file
    )