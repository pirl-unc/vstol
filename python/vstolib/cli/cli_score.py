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
and run 'score' command.
"""


import argparse
import numpy as np
import pandas as pd
import multiprocessing as mp
import logging
from ..constants import *
from ..default import *
from ..logging import get_logger
from ..main import score
from ..utilities import str2bool
from ..variants_list import VariantsList


logger = get_logger(__name__)


def add_cli_score_arg_parser(
        sub_parsers: argparse._SubParsersAction
) -> argparse._SubParsersAction:
    """
    Adds 'score' parser.
    """
    parser = sub_parsers.add_parser(
        'score',
        help='Calculates average alignment score for each breakpoint.'
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--tsv-file", '-i',
        dest="tsv_file",
        type=str,
        required=True,
        help="Variants list TSV file. "
             "This TSV file must follow VSTOL's TSV format for this command to work properly."
    )
    parser_required.add_argument(
        "--bam-file", '-b',
        dest="bam_file",
        type=str,
        required=True,
        help="BAM file."
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
        "--window",
        dest="window",
        type=int,
        required=False,
        default=WINDOW,
        help="Breakpoint window (default: %i). The average alignment score will be calculated "
             "with reads mapped upstream and downstream of this window from each breakpoint."
             % WINDOW
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
    parser_optional.add_argument(
        "--num-threads",
        dest="num_threads",
        type=int,
        default=NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)." % NUM_THREADS
    )

    parser.set_defaults(which='score')
    return sub_parsers


def run_cli_score_from_parsed_args(args: argparse.Namespace):
    """
    Run 'score' command using parameters from parsed arguments.

    Parameters:
        args    :   argparse.Namespace object with the following variables:
                    tsv_file
                    bam_file
                    output_tsv_file
                    window
                    gzip
                    num_threads
    """
    # Step 1. Load variants
    variants_list = VariantsList.read_tsv_file(tsv_file=args.tsv_file)

    # Step 2. Calculate average alignment score for each breakpoint
    logger.info("Started calculating average alignment score for each breakpoint.")
    variants_list = score(
        variants_list=variants_list,
        bam_file=args.bam_file,
        window=args.window,
        num_threads=args.num_threads
    )
    logger.info("Finished calculating average alignment score for each breakpoint.")

    # Step 3. Write to a TSV file
    df_variants = variants_list.to_dataframe()
    df_variants.sort_values(['variant_id'], inplace=True)
    if args.gzip:
        if args.output_tsv_file.endswith(".gz") == False:
            args.output_tsv_file = args.output_tsv_file + '.gz'
        df_variants.to_csv(args.output_tsv_file, sep='\t', index=False, compression='gzip')
    else:
        df_variants.to_csv(args.output_tsv_file, sep='\t', index=False)
