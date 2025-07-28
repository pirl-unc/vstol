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
and run 'diff' command.
"""


import argparse
import pandas as pd
from ..default import *
from ..logging import get_logger
from ..main import diff
from ..utilities import str2bool


logger = get_logger(__name__)


def add_cli_diff_arg_parser(
        sub_parsers: argparse._SubParsersAction
) -> argparse._SubParsersAction:
    """
    Adds 'diff' parser.
    """
    parser = sub_parsers.add_parser(
        'diff',
        help='Diff variants (outputs a variants list with variants specific to the target variants list.'
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--target-tsv-file", '-i',
        dest="target_tsv_file",
        required=True,
        help="Target variants list TSV file. "
             "This TSV file must follow the Occam's Variant Grammar format for this command to "
             "work properly."
    )
    parser_required.add_argument(
        "--query-tsv-file", '-q',
        dest="query_tsv_files",
        nargs='+',
        required=True,
        help="Query variants list TSV file. "
             "This TSV file must follow the Occam's Variant Grammar format for this command to "
             "work properly."
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
        "--num-processes",
        dest="num_processes",
        type=int,
        default=NUM_PROCESSES,
        required=False,
        help="Number of processes (default: %i)." % NUM_PROCESSES
    )
    parser_optional.add_argument(
        "--match-both-positions",
        dest="match_both_positions",
        type=str2bool,
        required=False,
        default=MATCH_BOTH_POSITIONS,
        help="If 'yes', two variants are considered an intersect "
             "if both pairs of breakpoints match ((position_1==position_1 AND position_2==position_2) OR (position_1==position_2 AND position_2==position_1)). "
             "If 'no', two variants are considered an intersect "
             "if one of the position pairs matches (position_1==position_1 OR position_1==position_2 OR position_2==position_1 OR position_2==position_2). Default: %s"
             % MATCH_BOTH_POSITIONS
    )
    parser_optional.add_argument(
        "--max-breakpoint-distance",
        dest="max_breakpoint_distance",
        type=int,
        required=False,
        default=MAX_BREAKPOINT_DISTANCE,
        help="Maximum breakpoint distance (default: %i). Two variants within "
             "this distance (inclusive) will be identified as the same variant."
             % MAX_BREAKPOINT_DISTANCE
    )
    parser_optional.add_argument(
        "--match-operation-types",
        dest="match_operation_types",
        type=str2bool,
        required=False,
        default=MATCH_OPERATION_TYPES,
        help="If 'yes', two variants are considered an intersect "
             "if the operation types are the same. Default: %s."
             % MATCH_OPERATION_TYPES
    )

    parser.set_defaults(which='diff')
    return sub_parsers


def run_cli_diff_from_parsed_args(args: argparse.Namespace):
    """
    Run 'diff' command using parameters from parsed arguments.

    Parameters:
        args    :   argparse.Namespace object with the following variables:
                    target_tsv_file
                    query_tsv_file
                    output_tsv_file
                    num_processes
                    match_both_positions
                    max_breakpoint_distance
                    match_operation_types
    """
    # Step 1. Load variants lists
    df_variants = pd.read_csv(args.target_tsv_file, sep='\t', low_memory=False, memory_map=True)
    df_variants_list = []
    for tsv_file in args.query_tsv_files:
        df_variants_list.append(pd.read_csv(tsv_file, sep='\t', low_memory=False, memory_map=True))

    # Step 2. Subtract variants lists
    df_filtered = diff(
        df_variants=df_variants,
        df_variants_list=df_variants_list,
        match_both_positions=args.match_both_positions,
        max_breakpoint_distance=args.max_breakpoint_distance,
        match_operation_types=args.match_operation_types,
        num_processes=args.num_processes
    )

    # Step 3. Write to a TSV file
    df_filtered.to_csv(args.output_tsv_file, sep='\t', index=False)
