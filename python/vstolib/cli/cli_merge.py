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
and run 'merge' command.
"""


import argparse
import numpy as np
import pandas as pd
import multiprocessing as mp
import logging
from ..constants import *
from ..default import *
from ..logging import get_logger
from ..main import merge
from ..utilities import str2bool
from ..variants_list import VariantsList


logger = get_logger(__name__)


def add_cli_merge_arg_parser(
        sub_parsers: argparse._SubParsersAction
) -> argparse._SubParsersAction:
    """
    Adds 'merge' parser.
    """
    parser = sub_parsers.add_parser(
        'merge',
        help='Merge variants (take union and deduplicate). '
             'Note that for two variant calls to be merged into'
             'one variant, (1) both breakpoints (position_1 and position_2) '
             'of the two calls must be nearby and the variant types '
             'must be of the same group. Note that INS=DUP and BND=INV=TRA for '
             'grouping purposes.'
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--tsv-file", '-i',
        dest="tsv_file",
        action='append',
        required=True,
        help="Variants list TSV file. "
             "This TSV file must follow VSTOL's TSV format for this command to work properly."
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
        "--num-threads", '-t',
        dest="num_threads",
        type=int,
        default=NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)." % NUM_THREADS
    )
    parser_optional.add_argument(
        "--max-neighbor-distance", '-d',
        dest="max_neighbor_distance",
        type=int,
        required=False,
        default=MERGE_MAX_NEIGHBOR_DISTANCE,
        help="Maximum neighbor distance (default: %i). Two variant calls within "
             "this distance (inclusive) will be merged into one variant."
             % MERGE_MAX_NEIGHBOR_DISTANCE
    )
    parser_optional.add_argument(
        "--match-all-breakpoints", '-b',
        dest="match_all_breakpoints",
        type=str2bool,
        required=False,
        default=MERGE_MATCH_ALL_BREAKPOINTS,
        help="If 'yes', two VariantCall objects are considered an intersect "
             "if both pairs of breakpoints match (start1==start2 AND end1==end2). "
             "If 'no', two VariantCall objects are considered an intersect "
             "if one of the breakpoint pairs matches (start1==start2 OR end1==end2). Default: %s"
             % MERGE_MATCH_ALL_BREAKPOINTS
    )
    parser_optional.add_argument(
        "--match-variant-types", '-m',
        dest="match_variant_types",
        type=str2bool,
        required=False,
        default=MERGE_MATCH_VARIANT_TYPES,
        help="If 'yes', two VariantCall objects are considered an intersect "
             "if the (super) variant types are the same. "
             "If 'no', two VariantCall objects are considered an intersect "
             "even if the (super) variant types are different (default: %s)."
             % MERGE_MATCH_VARIANT_TYPES
    )
    parser.set_defaults(which='merge')
    return sub_parsers


def run_cli_merge_from_parsed_args(args: argparse.Namespace):
    """
    Run 'merge' command using parameters from parsed arguments.

    Parameters:
        args    :   argparse.Namespace object with the following variables:
                    tsv_file
                    output_tsv_file
                    num_threads
                    max_neighbor_distance
                    match_all_breakpoints
                    match_all_breakpoints
    """
    # Step 1. Load variants lists
    logger.info("Started reading all TSV files")
    pool = mp.Pool(processes=args.num_threads)
    async_results = []
    for tsv_file in args.tsv_file:
        async_results.append(pool.apply_async(VariantsList.read_tsv_file, args=(tsv_file,True,False)))
    pool.close()
    pool.join()
    variants_lists = [async_result.get() for async_result in async_results]
    logger.info("Finished reading all TSV files")

    # Step 2. Merge variants lists
    logger.info("Started merging all variants into one list")
    variants_list = merge(
        variants_lists=variants_lists,
        num_threads=args.num_threads,
        max_neighbor_distance=args.max_neighbor_distance,
        match_all_breakpoints=args.match_all_breakpoints,
        match_variant_types=args.match_variant_types
    )
    logger.info("Finished merging all variants into one list")
    logger.info('%i variants and %i variant calls in merged variants list' %
                (len(variants_list.variant_ids), len(variants_list.variant_call_ids)))

    # Step 3. Write to a TSV file
    df_variants_list = variants_list.to_dataframe()
    df_variants_list.sort_values(['variant_id'], inplace=True)
    df_variants_list.to_csv(
        args.output_tsv_file,
        sep='\t',
        index=False
    )
