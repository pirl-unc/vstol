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
and run 'vcf2tsv' command.
"""


import argparse
from ..constants import *
from ..default import *
from ..logging import get_logger
from ..main import vcf2tsv
from ..utilities import str2bool
from ..vcf.common import read_vcf_file


logger = get_logger(__name__)


def add_cli_vcf2tsv_arg_parser(
        sub_parsers: argparse._SubParsersAction
) -> argparse._SubParsersAction:
    """
    Add 'vcf2tsv' parser.
    """
    parser = sub_parsers.add_parser('vcf2tsv', help='Convert a VCF file to a TSV file.')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--vcf-file", '-i',
        dest="vcf_file",
        type=str,
        required=True,
        help="Input VCF file."
    )
    parser_required.add_argument(
        "--variant-calling-method", '-m',
        dest="variant_calling_method",
        type=str,
        required=True,
        choices=VariantCallingMethods.ALL,
        help="Variant calling method. "
             "Allowed options: %s."
             % (', '.join(VariantCallingMethods.ALL))
    )
    parser_required.add_argument(
        "--sequencing-platform", '-p',
        dest="sequencing_platform",
        type=str,
        required=True,
    )
    parser_required.add_argument(
        "--source-id", '-s',
        dest="source_id",
        type=str,
        required=True,
        help="Source ID (e.g. patient ID or cell line ID)."
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
        "--case-id",
        dest="case_id",
        type=str,
        required=False,
        help="Case ID. This parameter must be specified if --variant-calling-method is strelka2-somatic. "
             "This is the tumor (case) ID that corresponds to 'TUMOR'."
    )
    parser_optional.add_argument(
        "--control-id",
        dest="control_id",
        type=str,
        required=False,
        help="Control ID. This parameter must be specified if --variant-calling-method is strelka2-somatic. "
             "This is the normal (control) ID that corresponds to 'NORMAL'."
    )
    parser_optional.add_argument(
        "--gzip",
        dest="gzip",
        type=str2bool,
        required=False,
        default=VCF2TSV_GZIP,
        help="If 'yes', gzip the output TSV file (default: %s)."
             % VCF2TSV_GZIP
    )

    parser.set_defaults(which='vcf2tsv')
    return sub_parsers


def run_cli_vcf2tsv_from_parsed_args(args: argparse.Namespace):
    """
    Run 'vcf2tsv' command using parameters from parsed arguments.

    Parameters:
        args    :   argparse.Namespace object with the following variables:
                    vcf_file
                    variant_calling_method
                    sequencing_platform
                    source_id
                    output_tsv_file
                    gzip
    """
    if args.variant_calling_method == VariantCallingMethods.STRELKA2_SOMATIC:
        if args.case_id is None or args.control_id is None:
            raise Exception("The parameters --case-id and --control-id must be "
                            "specified when --variant-calling-method strelka2-somatic")
    if args.variant_calling_method == VariantCallingMethods.SAVANA:
        if args.case_id is None or args.control_id is None:
            raise Exception("The parameters --case-id and --control-id must be "
                            "specified when --variant-calling-method savana")
    if args.variant_calling_method == VariantCallingMethods.SEVERUS:
        if args.case_id is None:
            raise Exception("The parameter --case-id must be "
                            "specified when --variant-calling-method severus")
    if args.variant_calling_method == VariantCallingMethods.SVISIONPRO:
        if args.case_id is None or args.control_id is None:
            raise Exception("The parameters --case-id and --control-id must be "
                            "specified when --variant-calling-method svisionpro")
    df_vcf = read_vcf_file(vcf_file=args.vcf_file)
    variants_list = vcf2tsv(
        df_vcf=df_vcf,
        source_id=args.source_id,
        variant_calling_method=args.variant_calling_method,
        sequencing_platform=args.sequencing_platform,
        case_id=args.case_id,
        control_id=args.control_id
    )

    df_variants = variants_list.to_dataframe()

    if args.gzip:
        if args.output_tsv_file.endswith(".gz") == False:
            args.output_tsv_file = args.output_tsv_file + '.gz'
        df_variants.to_csv(args.output_tsv_file, sep='\t', index=False, compression='gzip')
    else:
        df_variants.to_csv(args.output_tsv_file, sep='\t', index=False)
