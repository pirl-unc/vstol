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
import numpy as np
import pandas as pd
import multiprocessing as mp
import logging
from ..constants import *
from ..default import *
from ..logging import get_logger
from ..main import diff
from ..utilities import str2bool
from ..variants_list import VariantsList


logger = get_logger(__name__)


def add_cli_diff_arg_parser(
        sub_parsers: argparse._SubParsersAction
) -> argparse._SubParsersAction:
    """
    Adds 'diff' parser.
    """
    parser = sub_parsers.add_parser(
        'diff',
        help='Diff variants (outputs a variants list with variant calls specific to the target variants list. '
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
             "This TSV file must follow Velo's TSV format for this command to "
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
        default=DIFF_MAX_NEIGHBOR_DISTANCE,
        help="Maximum neighbor distance (default: %i). Two variant calls within "
             "this distance (inclusive) will be identified as the same variant."
             % DIFF_MAX_NEIGHBOR_DISTANCE
    )
    parser_optional.add_argument(
        "--match-all-breakpoints",
        dest="match_all_breakpoints",
        type=str2bool,
        required=False,
        default=DIFF_MATCH_ALL_BREAKPOINTS,
        help="If 'yes', two VariantCall objects are considered an intersect "
             "if both pairs of breakpoints match (start1==start2 AND end1==end2). "
             "If 'no', two VariantCall objects are considered an intersect "
             "if one of the breakpoint pairs matches (start1==start2 OR end1==end2). Default: %s"
             % DIFF_MATCH_ALL_BREAKPOINTS
    )
    parser_optional.add_argument(
        "--match-variant-types",
        dest="match_variant_types",
        type=str2bool,
        required=False,
        default=DIFF_MATCH_VARIANT_TYPES,
        help="If 'yes', two VariantCall objects are considered an intersect "
             "if the (super) variant types are the same. "
             "If 'no', two VariantCall objects are considered an intersect "
             "even if the (super) variant types are different (default: %s)."
             % DIFF_MATCH_VARIANT_TYPES
    )
    parser_required.add_argument(
        "--query-bam-file",
        dest="query_bam_file",
        action='append',
        required=False,
        help="Query BAM file. If specified, presence of target variants "
             "is checked in these BAM files (recommended for somatic or "
             "case-specific variant identification). "
             "Every BAM file must have CS and MD tags for this "
             "command to work properly."
    )
    parser_optional.add_argument(
        "--gzip",
        dest="gzip",
        type=str2bool,
        required=False,
        default=DIFF_GZIP,
        help="If 'yes', gzip the output TSV file (default: %s)."
             % DIFF_GZIP
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
                    num_threads
                    max_neighbor_distance
                    match_all_breakpoints
                    match_variant_types
                    query_bam_file
                    gzip
    """
    # Step 1. Load the target variants list
    logger.info("Started reading the target variants list TSV file")
    target_variants_list = VariantsList.read_tsv_file(tsv_file=args.target_tsv_file)
    logger.info("Finished reading the target variants list TSV file")

    # Step 2. Diff each query variants list one at a time
    logger.info("Started diffing")
    for tsv_file in args.query_tsv_files:
        logger.info("Loading %s" % tsv_file)
        query_variants_list = VariantsList.read_tsv_file(tsv_file=tsv_file)
        target_variants_list = diff(
            target_variants_list=target_variants_list,
            query_variants_lists=[query_variants_list],
            num_threads=args.num_threads,
            max_neighbor_distance=args.max_neighbor_distance,
            match_all_breakpoints=args.match_all_breakpoints,
            match_variant_types=args.match_variant_types
        )
    logger.info("Finished diffing")

    # Step 3. Write to a TSV file
    df_variants = target_variants_list.to_dataframe()
    df_variants.sort_values(['variant_id'], inplace=True)
    if args.gzip:
        if args.output_tsv_file.endswith(".gz") == False:
            args.output_tsv_file = args.output_tsv_file + '.gz'
        df_variants.to_csv(args.output_tsv_file, sep='\t', index=False, compression='gzip')
    else:
        df_variants.to_csv(args.output_tsv_file, sep='\t', index=False)

