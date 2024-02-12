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
import numpy as np
import pandas as pd
import multiprocessing as mp
import logging
from ..constants import *
from ..default import *
from ..genomic_ranges_list import GenomicRangesList
from ..logging import get_logger
from ..main import overlap
from ..variants_list import VariantsList


logger = get_logger(__name__)


def add_cli_overlap_arg_parser(
        sub_parsers: argparse._SubParsersAction
) -> argparse._SubParsersAction:
    """
    Adds 'overlap' parser.
    """
    parser = sub_parsers.add_parser(
        'overlap',
        help='Overlap variants in a list of genomic ranges.'
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
        "--num-threads",
        dest="num_threads",
        type=int,
        default=NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)." % NUM_THREADS
    )
    parser_optional.add_argument(
        "--padding",
        dest="padding",
        type=int,
        required=False,
        default=OVERLAP_PADDING,
        help="Padding to apply to each breakpoint (default: %i)."
             % OVERLAP_PADDING
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
                    num_threads
                    padding
    """
    # Step 1. Load variants list
    logger.info("Loading variants list")
    variants_list = VariantsList.read_tsv_file(tsv_file=args.tsv_file)

    # Step 2. Load regions list
    logger.info("Loading genomic ranges list")
    genomic_ranges_list = GenomicRangesList.read_tsv_file(
        tsv_file=args.regions_tsv_file
    )

    # Step 3. Identify overlaps
    logger.info("Started identifying overlaps")
    variants_list = overlap(
        variants_list=variants_list,
        genomic_ranges_list=genomic_ranges_list,
        padding=args.padding,
        num_threads=args.num_threads
    )
    logger.info("Finished identifying overlaps")
    logger.info('%i variants and %i variant calls overlap' %
                (len(variants_list.variant_ids), len(variants_list.variant_call_ids)))

    # Step 3. Write to a TSV file
    df_variants_list = variants_list.to_dataframe()
    df_variants_list.sort_values(['variant_id'], inplace=True)
    df_variants_list.to_csv(
        args.output_tsv_file,
        sep='\t',
        index=False
    )
