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
The purpose of this python3 script is to implement the primary VSTOL command.
"""


import argparse
import vstolib
from typing import Tuple
from .cli_annotate import *
from .cli_collapse import *
from .cli_diff import *
from .cli_filter import *
from .cli_intersect import *
from .cli_merge import *
from .cli_overlap import *
from .cli_tsv2vcf import *
from .cli_vcf2tsv import *
from ..logging import get_logger


logger = get_logger(__name__)


def init_arg_parser() -> Tuple[argparse.ArgumentParser, argparse._SubParsersAction]:
    """
    Initialize the input argument parser.

    Returns:
        Tuple[argparse.ArgumentParser,argparse.ArgumentParser subparsers]
    """
    arg_parser = argparse.ArgumentParser(
        description="VSTOL: Variant Selection, Tabulation, and Operations Library."
    )
    arg_parser.add_argument(
        '--version', '-v',
        action='version',
        version='%(prog)s version ' + str(vstolib.__version__)
    )
    sub_parsers = arg_parser.add_subparsers(help='VSTOL sub-commands.')
    return arg_parser, sub_parsers


def run():
    # Step 1. Initialize argument parser
    arg_parser, sub_parsers = init_arg_parser()
    sub_parsers = add_cli_annotate_arg_parser(sub_parsers=sub_parsers)     # annotate
    sub_parsers = add_cli_collapse_arg_parser(sub_parsers=sub_parsers)     # collapse
    sub_parsers = add_cli_diff_arg_parser(sub_parsers=sub_parsers)         # diff
    sub_parsers = add_cli_filter_arg_parser(sub_parsers=sub_parsers)       # filter
    sub_parsers = add_cli_intersect_arg_parser(sub_parsers=sub_parsers)    # intersect
    sub_parsers = add_cli_merge_arg_parser(sub_parsers=sub_parsers)        # merge
    sub_parsers = add_cli_overlap_arg_parser(sub_parsers=sub_parsers)      # overlap
    sub_parsers = add_cli_tsv2vcf_arg_parser(sub_parsers=sub_parsers)      # tsv2vcf
    sub_parsers = add_cli_vcf2tsv_arg_parser(sub_parsers=sub_parsers)      # vcf2tsv
    args = arg_parser.parse_args()

    # Step 2. Execute function based on CLI arguments
    if args.which == 'annotate':
        run_cli_annotate_from_parsed_args(args=args)
    elif args.which == 'collapse':
        run_cli_collapse_from_parsed_args(args=args)
    elif args.which == 'diff':
        run_cli_diff_from_parsed_args(args=args)
    elif args.which == 'filter':
        run_cli_filter_from_parsed_args(args=args)
    elif args.which == 'intersect':
        run_cli_intersect_from_parsed_args(args=args)
    elif args.which == 'merge':
        run_cli_merge_from_parsed_args(args=args)
    elif args.which == 'overlap':
        run_cli_overlap_from_parsed_args(args=args)
    elif args.which == 'tsv2vcf':
        run_cli_tsv2vcf_from_parsed_args(args=args)
    elif args.which == 'vcf2tsv':
        run_cli_vcf2tsv_from_parsed_args(args=args)
    else:
        raise Exception("Invalid command: %s" % args.which)
