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
and run 'overlap' command.
"""


import argparse
import pandas as pd
from ..default import *
from ..logging import get_logger
from ..main import overlap


logger = get_logger(__name__)


def add_cli_overlap_arg_parser(
        sub_parsers: argparse._SubParsersAction
) -> argparse._SubParsersAction:
    """
    Adds 'overlap' parser.
    """
    parser = sub_parsers.add_parser(
        'overlap',
        help='Identify a set of variants overlapping genomic ranges.'
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--tsv-file", '-i',
        dest="tsv_file",
        type=str,
        required=True,
        help="Variant list TSV file. Expected columns: 'id', 'chromosome_1', position_1', 'operation_1', 'chromosome_2', 'position_2', 'operation_2'."
    )
    parser_required.add_argument(
        "--regions-tsv-file", '-r',
        dest="regions_tsv_file",
        type=str,
        required=False,
        help="TSV file of genomic regions. "
             "Variants with any variant calls where breakpoints are near or within the regions "
             "in this file will be considered an overlap. "
             "Expected headers: 'chromosome', 'start', 'end'."
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
        "--buffer",
        dest="buffer",
        type=int,
        default=BUFFER,
        required=False,
        help="Buffer (default: %i)." % BUFFER
    )

    parser.set_defaults(which='overlap')
    return sub_parsers


def run_cli_overlap_from_parsed_args(args: argparse.Namespace):
    """
    Run 'overlap' command using parameters from parsed arguments.

    Parameters:
        args    :   argparse.Namespace object with the following variables:
                    tsv_file
                    ranges_tsv_file
                    output_tsv_file
                    buffer
    """
    df_variants = pd.read_csv(args.tsv_file, sep='\t', low_memory=False, memory_map=True)
    df_genomic_ranges = pd.read_csv(args.regions_tsv_file, sep='\t', low_memory=False, memory_map=True)
    df_variants_overlapping = overlap(
        df_variants=df_variants,
        df_genomic_ranges=df_genomic_ranges,
        buffer=args.buffer
    )
    df_variants_overlapping.to_csv(args.output_tsv_file, sep='\t', index=False)
