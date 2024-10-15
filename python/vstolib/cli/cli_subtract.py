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
and run 'subtract' command.
"""


import argparse
import numpy as np
import pandas as pd
import multiprocessing as mp
import logging
from ..constants import *
from ..default import *
from ..logging import get_logger
from ..main import subtract
from ..utilities import str2bool
from ..variants_list import VariantsList


logger = get_logger(__name__)


def add_cli_subtract_arg_parser(
        sub_parsers: argparse._SubParsersAction
) -> argparse._SubParsersAction:
    """
    Adds 'subtract' parser.
    """
    parser = sub_parsers.add_parser(
        'subtract',
        help='Subtract variants (outputs a variants list with variant calls specific to the target variants list. '
             'Note that for two variant calls to be considered the same variant, '
             '(1) both breakpoints (position_1 and position_2) '
             'of the two calls must be nearby and the variant types '
             'must be of the same group. Note that INS=DUP and BND=INV=TRA for '
             'grouping purposes.'
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--target-tsv-file", '-i',
        dest="target_tsv_file",
        required=True,
        help="Target variants list TSV file. "
             "This TSV file must follow VSTOL's TSV format for this command to "
             "work properly."
    )
    parser_required.add_argument(
        "--query-tsv-file", '-q',
        dest="query_tsv_files",
        action='append',
        required=True,
        help="Query variants list TSV file. "
             "This TSV file must follow VSTOL's TSV format for this command to "
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
        "--num-threads",
        dest="num_threads",
        type=int,
        default=NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)." % NUM_THREADS
    )
    parser_optional.add_argument(
        "--max-neighbor-distance",
        dest="max_neighbor_distance",
        type=int,
        required=False,
        default=MAX_NEIGHBOR_DISTANCE,
        help="Maximum neighbor distance (default: %i). Two variant calls within "
             "this distance (inclusive) will be identified as the same variant."
             % MAX_NEIGHBOR_DISTANCE
    )
    parser_optional.add_argument(
        "--match-all-breakpoints",
        dest="match_all_breakpoints",
        type=str2bool,
        required=False,
        default=MATCH_ALL_BREAKPOINTS,
        help="If 'yes', two VariantCall objects are considered an intersect "
             "if both pairs of breakpoints match (start1==start2 AND end1==end2). "
             "If 'no', two VariantCall objects are considered an intersect "
             "if one of the breakpoint pairs matches (start1==start2 OR end1==end2). Default: %s"
             % MATCH_ALL_BREAKPOINTS
    )
    parser_optional.add_argument(
        "--match-variant-types",
        dest="match_variant_types",
        type=str2bool,
        required=False,
        default=MATCH_VARIANT_TYPES,
        help="If 'yes', two VariantCall objects are considered an intersect "
             "if the (super) variant types are the same. "
             "If 'no', two VariantCall objects are considered an intersect "
             "even if the (super) variant types are different (default: %s)."
             % MATCH_VARIANT_TYPES
    )
    parser_optional.add_argument(
        "--min-ins-size-overlap",
        dest="min_ins_size_overlap",
        type=int,
        required=False,
        default=MIN_INS_SIZE_OVERLAP,
        help="Minimum insertion size overlap (default: %f)."
             % MIN_INS_SIZE_OVERLAP
    )
    parser_optional.add_argument(
        "--min-del-size-overlap",
        dest="min_del_size_overlap",
        type=int,
        required=False,
        default=MIN_DEL_SIZE_OVERLAP,
        help="Minimum deletion size overlap (default: %f)."
             % MIN_DEL_SIZE_OVERLAP
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

    parser.set_defaults(which='subtract')
    return sub_parsers


def run_cli_subtract_from_parsed_args(args: argparse.Namespace):
    """
    Run 'subtract' command using parameters from parsed arguments.

    Parameters:
        args    :   argparse.Namespace object with the following variables:
                    target_tsv_file
                    query_tsv_file
                    output_tsv_file
                    num_threads
                    max_neighbor_distance
                    match_all_breakpoints
                    match_variant_types
                    min_ins_size_overlap
                    min_del_size_overlap
                    gzip
    """
    # Step 1. Load the target variants list
    logger.info("Started reading the target variants list TSV file")
    target_variants_list = VariantsList.read_tsv_file(tsv_file=args.target_tsv_file)
    logger.info("Finished reading the target variants list TSV file")

    # Step 2. Diff each query variants list one at a time
    logger.info("Started subtracting")
    for tsv_file in args.query_tsv_files:
        logger.info("Loading %s" % tsv_file)
        query_variants_list = VariantsList.read_tsv_file(tsv_file=tsv_file)
        target_variants_list = subtract(
            target_variants_list=target_variants_list,
            query_variants_lists=[query_variants_list],
            num_threads=args.num_threads,
            max_neighbor_distance=args.max_neighbor_distance,
            match_all_breakpoints=args.match_all_breakpoints,
            match_variant_types=args.match_variant_types,
            min_ins_size_overlap=args.min_ins_size_overlap,
            min_del_size_overlap=args.min_del_size_overlap
        )
    logger.info("Finished subtracting")

    # Step 3. Write to a TSV file
    df_variants = target_variants_list.to_dataframe()
    df_variants.sort_values(['variant_id'], inplace=True)
    if args.gzip:
        if args.output_tsv_file.endswith(".gz") == False:
            args.output_tsv_file = args.output_tsv_file + '.gz'
        df_variants.to_csv(args.output_tsv_file, sep='\t', index=False, compression='gzip')
    else:
        df_variants.to_csv(args.output_tsv_file, sep='\t', index=False)

