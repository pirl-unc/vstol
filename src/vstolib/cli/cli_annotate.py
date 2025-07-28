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
and run 'annotate' command.
"""


import argparse
import pandas as pd
from ..constants import *
from ..default import *
from ..gencode import Gencode
from ..main import annotate


def add_cli_annotate_arg_parser(
        sub_parsers: argparse._SubParsersAction
) -> argparse._SubParsersAction:
    """
    Adds 'annotate' parser.
    """
    parser = sub_parsers.add_parser(
        'annotate',
        help='Annotate variants.'
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--tsv-file", '-i',
        dest="tsv_file",
        type=str,
        required=True,
        help="Input variants TSV file. Expected headers: "
             "'chromosome_1', 'position_1', "
             "'chromosome_2', 'position_2'"
    )
    parser_required.add_argument(
        "--annotator", '-a',
        dest="annotator",
        type=str,
        choices=Annotator.ALL,
        required=True,
        help="Annotator. Allowed options: %s."
             % (', '.join(f"'{item}'" for item in Annotator.ALL))
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
        help="Number of processes (default: %i)."
             % NUM_PROCESSES
    )
    parser_optional.add_argument(
        "--gencode-gtf-file",
        dest="gencode_gtf_file",
        type=str,
        required=False,
        help="GENCODE GTF file. "
             "This parameter must be supplied if "
             "--annotator is '%s'."
             % Annotator.GENCODE
    )
    parser_optional.add_argument(
        "--gencode-levels",
        dest="gencode_levels",
        nargs='+',
        type=int,
        required=False,
        default=[1],
        help="GENCODE gene levels (default: 1)."
             "This parameter must be supplied if "
             "--annotator is '%s'."
             % Annotator.GENCODE
    )
    parser_optional.add_argument(
        "--gencode-types",
        dest="gencode_types",
        nargs='+',
        type=str,
        required=False,
        default=[],
        help="GENCODE gene types (default: [])."
             "This parameter must be supplied if "
             "--annotator is '%s'."
             % Annotator.GENCODE
    )
    parser_optional.add_argument(
        "--gencode-version",
        dest="gencode_version",
        type=str,
        required=False,
        help="GENCODE version (e.g. 'v41'). "
             "This parameter must be supplied if "
             "--annotator is '%s'."
             % Annotator.GENCODE
    )
    parser_optional.add_argument(
        "--gencode-species",
        dest="gencode_species",
        type=str,
        required=False,
        help="GENCODE species (e.g. 'human'). "
             "This parameter must be supplied if "
             "--annotator is '%s'."
             % Annotator.GENCODE
    )

    parser.set_defaults(which='annotate')
    return sub_parsers


def run_cli_annotate_from_parsed_args(args: argparse.Namespace):
    """
    Run 'annotate' command using parameters from parsed arguments.

    Parameters:
        args    :   argparse.Namespace object with the following variables:
                    tsv_file
                    annotator
                    output_tsv_file
                    gencode_gtf_file
                    gencode_levels
                    gencode_types
                    gencode_version
                    gencode_species
    """
    if args.annotator == Annotator.GENCODE:
        assert args.gencode_gtf_file is not None
        assert args.gencode_version is not None
        assert args.gencode_species is not None
        assert args.gencode_types is not None
        assert args.gencode_version is not None
        annotator = Gencode(
            gtf_file=args.gencode_gtf_file,
            version=args.gencode_version,
            species=args.gencode_species,
            types=args.gencode_types,
            levels=args.gencode_levels
        )
    else:
        raise Exception('annotator must be one of %s' % Annotator)

    df_variants = pd.read_csv(args.tsv_file, sep='\t', low_memory=False, memory_map=True)
    df_variants_annotated = annotate(
        df_variants=df_variants,
        annotator=annotator,
        num_processes=args.num_processes
    )

    df_variants_annotated.to_csv(args.output_tsv_file, sep='\t', index=False)
