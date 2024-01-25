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
        default=DIFF_MAX_NEIGHBOR_DISTANCE,
        help="Maximum neighbor distance (default: %i). Two variant calls within "
             "this distance (inclusive) will be identified as the same variant."
             % DIFF_MAX_NEIGHBOR_DISTANCE
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
    """
    # Step 1. Load variants lists
    logger.info("Started reading all TSV files")
    target_variants_list = VariantsList.read_tsv_file(tsv_file=args.target_tsv_file)
    pool = mp.Pool(processes=args.num_threads)
    async_results = []
    for tsv_file in args.query_tsv_files:
        async_results.append(pool.apply_async(VariantsList.read_tsv_file, args=(tsv_file,True,False)))
    pool.close()
    pool.join()
    query_variants_lists = [async_result.get() for async_result in async_results]
    logger.info("Finished reading all TSV files")

    # Step 2. Diff variants lists
    logger.info("Started diffing variants")
    variants_list = diff(
        target_variants_list=target_variants_list,
        query_variants_lists=query_variants_lists,
        num_threads=args.num_threads,
        max_neighbor_distance=args.max_neighbor_distance
    )
    logger.info("Finished diffing variants")

    # Step 3. Write to a TSV file
    df_variants_list = variants_list.to_dataframe()
    df_variants_list.sort_values(['variant_id'], inplace=True)
    df_variants_list.to_csv(
        args.output_tsv_file,
        sep='\t',
        index=False
    )
