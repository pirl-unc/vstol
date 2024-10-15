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
and run 'collapse' command.
"""


import argparse
from ..constants import *
from ..default import *
from ..logging import get_logger
from ..main import collapse
from ..utilities import str2bool
from ..variants_list import VariantsList


logger = get_logger(__name__)


def add_cli_collapse_arg_parser(
        sub_parsers: argparse._SubParsersAction
) -> argparse._SubParsersAction:
    """
    Add 'collapse' parser.
    """
    parser = sub_parsers.add_parser(
        'collapse',
        help='Collapses a VSTOL TSV file (with multiple VariantCall rows per Variant) '
             'into a TSV file with one representative VariantCall row per Variant '
             'to a TSV file.')
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
        "--sample-id", '-s',
        dest="sample_id",
        type=str,
        required=True,
        help="Sample ID to retain."
    )
    parser_required.add_argument(
        "--output-tsv-file", '-o',
        dest="output_tsv_file",
        type=str,
        required=True,
        help="Output TSV file."
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

    parser.set_defaults(which='collapse')
    return sub_parsers


def run_cli_collapse_from_parsed_args(args: argparse.Namespace):
    """
    Run 'collapse' command using parameters from parsed arguments.

    Parameters:
        args    :   argparse.Namespace object with the following variables:
                    tsv_file
                    sample_id
                    output_tsv_file
                    strategy
    """
    variants_list = VariantsList.read_tsv_file(tsv_file=args.tsv_file)
    variants_list = collapse(
        variants_list=variants_list,
        sample_id=args.sample_id,
        strategy=args.strategy
    )
    df_variants = variants_list.to_dataframe()
    df_variants.sort_values(['variant_id'], inplace=True)
    if args.gzip:
        if args.output_tsv_file.endswith(".gz") == False:
            args.output_tsv_file = args.output_tsv_file + '.gz'
        df_variants.to_csv(args.output_tsv_file, sep='\t', index=False, compression='gzip')
    else:
        df_variants.to_csv(args.output_tsv_file, sep='\t', index=False)
