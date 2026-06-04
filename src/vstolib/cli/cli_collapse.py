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
import pandas as pd
from ..default import *
from ..logging import get_logger
from ..main import collapse


logger = get_logger(__name__)


def add_cli_collapse_arg_parser(
        sub_parsers: argparse._SubParsersAction
) -> argparse._SubParsersAction:
    """
    Adds 'collapse' parser.
    """
    parser = sub_parsers.add_parser(
        'collapse',
        help='Collapse a merged variants TSV into one row per group_id, '
             'applying somatic filters.'
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--tsv-file",
        dest="tsv_file",
        type=str,
        required=True,
        help="Merged variants list TSV file (output of `vstol merge`). "
             "Must contain a 'group_id' column."
    )
    parser_required.add_argument(
        "--output-tsv-file",
        dest="output_tsv_file",
        type=str,
        required=True,
        help="Output TSV file."
    )
    parser_required.add_argument(
        "--tumor-sample-ids",
        dest="tumor_sample_ids",
        nargs='+',
        required=True,
        help="One or more tumor sample IDs. Used to pick the canonical "
             "variant representation (positions, read counts, AF) per group."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--normal-sample-ids",
        dest="normal_sample_ids",
        nargs='+',
        required=False,
        default=[],
        help="Zero or more normal sample IDs. If provided, any group "
             "containing a normal-sample row with > 0 alternate-allele "
             "reads at either breakpoint is dropped (normal contamination "
             "filter). Omit this argument when the merged TSV contains no "
             "normal sample (e.g. tumor-only callsets) - the contamination "
             "filter will be skipped."
    )
    parser_optional.add_argument(
        "--min-methods",
        dest="min_methods",
        type=int,
        required=False,
        default=MIN_METHODS,
        help="Minimum number of unique calling methods required per group "
             "(default: %i). Groups with fewer unique methods are dropped."
             % MIN_METHODS
    )

    parser.set_defaults(which='collapse')
    return sub_parsers


def run_cli_collapse_from_parsed_args(args: argparse.Namespace):
    # Step 1. Load merged variants list
    df_variants = pd.read_csv(args.tsv_file, sep='\t', low_memory=False, memory_map=True)

    # Step 2. Collapse
    df_collapsed = collapse(
        df_variants=df_variants,
        tumor_sample_ids=args.tumor_sample_ids,
        normal_sample_ids=args.normal_sample_ids,
        min_methods=args.min_methods
    )

    # Step 3. Write to a TSV file
    df_collapsed.to_csv(args.output_tsv_file, sep='\t', index=False)
