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
and run 'compare' command.
"""


import argparse
import multiprocessing as mp
from ..constants import *
from ..default import *
from ..logging import get_logger
from ..main import merge
from ..utilities import str2bool
from ..variants_list import VariantsList


logger = get_logger(__name__)


def add_cli_compare_arg_parser(
        sub_parsers: argparse._SubParsersAction
) -> argparse._SubParsersAction:
    """
    Adds 'compare' parser.
    """
    parser = sub_parsers.add_parser('compare', help='Compare variant lists.')
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
             "this distance (inclusive) will be merged into one variant."
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

    parser.set_defaults(which='compare')
    return sub_parsers


def run_cli_compare_from_parsed_args(args: argparse.Namespace):
    """
    Run 'compare' command using parameters from parsed arguments.

    Parameters:
        args    :   argparse.Namespace object with the following variables:
                    tsv_file
                    output_tsv_file
                    num_threads
                    max_neighbor_distance
                    match_all_breakpoints
                    match_variant_types
                    min_ins_size_overlap
                    min_del_size_overlap
                    gzip
    """
    assert(len(args.tsv_file), 2)

    # Step 1. Load variants lists
    logger.info("Started reading all TSV files")
    pool = mp.Pool(processes=args.num_threads)
    async_results = []
    for tsv_file in args.tsv_file:
        async_results.append(pool.apply_async(VariantsList.read_tsv_file, args=(tsv_file,False,True)))
    pool.close()
    pool.join()
    variants_lists = [async_result.get() for async_result in async_results]
    logger.info("Finished reading all TSV files")

    # Step 2. Merge variants lists
    logger.info("Started comparing all variants into one list")
    variants_list = merge(
        variants_lists=variants_lists,
        num_threads=args.num_threads,
        max_neighbor_distance=args.max_neighbor_distance,
        match_all_breakpoints=args.match_all_breakpoints,
        match_variant_types=args.match_variant_types,
        min_ins_size_overlap=args.min_ins_size_overlap,
        min_del_size_overlap=args.min_del_size_overlap
    )
    logger.info("Finished merging all variants into one list")
    logger.info('%i variants and %i variant calls in merged variants list' %
                (len(variants_list.variant_ids), len(variants_list.variant_call_ids)))

    # Step 3. Write to a TSV file
    df_variants = variants_list.to_dataframe()
    df_variants.sort_values(['variant_id'], inplace=True)
    if args.gzip:
        if args.output_tsv_file.endswith(".gz") == False:
            args.output_tsv_file = args.output_tsv_file + '.gz'
        df_variants.to_csv(args.output_tsv_file, sep='\t', index=False, compression='gzip')
    else:
        df_variants.to_csv(args.output_tsv_file, sep='\t', index=False)

