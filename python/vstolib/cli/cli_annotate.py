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
from ..constants import *
from ..default import *
from ..ensembl import Ensembl
from ..gencode import Gencode
from ..refseq import RefSeq
from ..main import annotate
from ..variants_list import VariantsList


def add_cli_annotate_arg_parser(
        sub_parsers: argparse._SubParsersAction
) -> argparse._SubParsersAction:
    """
    Adds 'annotate' parser.
    """
    parser = sub_parsers.add_parser(
        'annotate',
        help='Annotate a variants list.'
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
             "'variant_id', 'variant_call_id', "
             "'chromosome_1', 'position_1', "
             "'chromosome_2', 'position_2', "
             "'reference_allele', 'variant_allele', "
             "'variant_type'"
    )
    parser_required.add_argument(
        "--annotator", '-a',
        dest="annotator",
        type=str,
        choices=Annotators.ALL,
        required=True,
        help="Annotator. Allowed options: %s."
             % (', '.join(f"'{item}'" for item in Annotators.ALL))
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
        help="Number of threads (default: %i)."
             % NUM_THREADS
    )
    parser_optional.add_argument(
        "--ensembl-release",
        dest="ensembl_release",
        type=int,
        required=False,
        help="Ensembl release version number "
             "(e.g. 75 for GRCh37 or 106 for GRCh38). "
             "This parameter must be supplied if "
             "--annotator is '%s'. "
             "Please make sure the specified "
             "ensembl version is installed using pyensembl."
             % Annotators.ENSEMBL
    )
    parser_optional.add_argument(
        "--ensembl-species",
        dest="ensembl_species",
        type=str,
        required=False,
        help="Ensembl species "
             "(e.g. 'human' or 'mouse'). "
             "This parameter must be supplied if "
             "--annotator is '%s'. "
             "Please make sure the specified "
             "ensembl species is installed using pyensembl."
             % Annotators.ENSEMBL
    )
    parser_optional.add_argument(
        "--gencode-gtf-file",
        dest="gencode_gtf_file",
        type=str,
        required=False,
        help="GENCODE GTF file. "
             "This parameter must be supplied if "
             "--annotator is '%s'."
             % Annotators.GENCODE
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
             % Annotators.GENCODE
    )
    parser_optional.add_argument(
        "--gencode-types",
        dest="gencode_types",
        nargs='+',
        type=str,
        required=False,
        default=['protein_coding'],
        help="GENCODE gene types (default: protein_coding)."
             "This parameter must be supplied if "
             "--annotator is '%s'."
             % Annotators.GENCODE
    )
    parser_optional.add_argument(
        "--gencode-version",
        dest="gencode_version",
        type=str,
        required=False,
        help="GENCODE version (e.g. 'v41'). "
             "This parameter must be supplied if "
             "--annotator is '%s'."
             % Annotators.GENCODE
    )
    parser_optional.add_argument(
        "--gencode-species",
        dest="gencode_species",
        type=str,
        required=False,
        help="GENCODE species (e.g. 'human'). "
             "This parameter must be supplied if "
             "--annotator is '%s'."
             % Annotators.GENCODE
    )
    parser_optional.add_argument(
        "--refseq-gtf-file",
        dest="refseq_gtf_file",
        type=str,
        required=False,
        help="RefSeq GTF file. "
             "This parameter must be supplied if "
             "--annotator is '%s'."
             % Annotators.REFSEQ
    )
    parser_optional.add_argument(
        "--refseq-assembly-report-txt-file",
        dest="refseq_assembly_report_txt_file",
        type=str,
        required=False,
        help="RefSeq assembly report TXT file. "
             "This parameter must be supplied if "
             "--annotator is '%s'."
             % Annotators.REFSEQ
    )
    parser_optional.add_argument(
        "--refseq-version",
        dest="refseq_version",
        type=str,
        required=False,
        help="RefSeq version (e.g. 'v110'). "
             "This parameter must be supplied if "
             "--annotator is '%s'."
             % Annotators.REFSEQ
    )
    parser_optional.add_argument(
        "--refseq-species",
        dest="refseq_species",
        type=str,
        required=False,
        help="RefSeq species (e.g. 'human'). "
             "This parameter must be supplied if "
             "--annotator is '%s'."
             % Annotators.REFSEQ
    )

    # parser_optional.add_argument(
    #     "--annovar-path", '-ar',
    #     dest="annovar_path",
    #     type=str,
    #     required=False,
    #     help="ANNOVAR directory path. "
    #          "This parameter must be supplied if "
    #          "--annotator is '%s'."
    #          % Annotators.ANNOVAR
    # )
    # parser_optional.add_argument(
    #     "--annovar-perl-path", '-ap',
    #     dest="annovar_perl_path",
    #     type=str,
    #     required=False,
    #     help="Perl path. "
    #          "This parameter must be supplied if "
    #          "--annotator is '%s'."
    #          % Annotators.ANNOVAR
    # )
    # parser_optional.add_argument(
    #     "--annovar-db-path", '-ad',
    #     dest="annovar_db_path",
    #     type=str,
    #     required=False,
    #     help="ANNOVAR <human>db/ directory path. "
    #          "This parameter must be supplied if "
    #          "--annotator is '%s'."
    #          % Annotators.ANNOVAR
    # )
    # parser_optional.add_argument(
    #     "--annovar-output-avinput-file", '-af',
    #     dest="annovar_output_avinput_file",
    #     type=str,
    #     required=False,
    #     help="Output ANNOVAR .avinput file. "
    #          "This parameter must be supplied if "
    #          "--annotator is '%s'."
    #          % Annotators.ANNOVAR
    # )
    # parser_optional.add_argument(
    #     "--annovar-genome-assembly", '-aa',
    #     dest="annovar_genome_assembly",
    #     type=str,
    #     required=False,
    #     help="ANNOVAR genome assembly. "
    #          "This parameter must be supplied if "
    #          "--annotator is '%s'."
    #          % Annotators.ANNOVAR
    # )
    # parser_optional.add_argument(
    #     "--annovar-protocol", '-at',
    #     dest="annovar_protocol",
    #     type=str,
    #     required=False,
    #     default=','.join(list(ANNOVAR_PROTOCOL_OPERATION.keys())),
    #     help="ANNOVAR protocol (e.g. 'refGene,exac03'). "
    #          "This parameter must be supplied if "
    #          "--annotator is '%s' (default: '%s')."
    #          % (Annotators.ANNOVAR, ','.join(list(ANNOVAR_PROTOCOL_OPERATION.keys())))
    # )
    # parser_optional.add_argument(
    #     "--annovar-operation", '-ao',
    #     dest="annovar_operation",
    #     type=str,
    #     required=False,
    #     default=','.join(list(ANNOVAR_PROTOCOL_OPERATION.values())),
    #     help="ANNOVAR protocol (e.g. 'g,f'). "
    #          "This parameter must be supplied if "
    #          "--annotator is '%s' (default: '%s')."
    #          % (Annotators.ANNOVAR, ','.join(list(ANNOVAR_PROTOCOL_OPERATION.values())))
    # )
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
                    ensembl_release
                    ensembl_species
                    gencode_gtf_file
                    gencode_levels
                    gencode_types
                    gencode_version
                    gencode_species
                    refseq_gtf_file
                    refseq_assembly_report_txt_file
                    refseq_version
                    refseq_species
                    annovar_path
                    annovar_perl_path
                    annovar_humandb_path
                    annovar_output_avinput_file
                    annovar_genome_assembly
                    annovar_protocol
                    annovar_operation
    """
    # Step 1. Load variants
    variants_list = VariantsList.read_tsv_file(tsv_file=args.tsv_file)

    # Step 2. Load annotation data
    if args.annotator == Annotators.ENSEMBL:
        annotator = Ensembl(
            release=args.ensembl_release,
            species=args.ensembl_species
        )
    elif args.annotator == Annotators.GENCODE:
        annotator = Gencode(
            gtf_file=args.gencode_gtf_file,
            version=args.gencode_version,
            species=args.gencode_species,
            types=args.gencode_types,
            levels=args.gencode_levels
        )
    else:
        raise Exception('Unknown annotation source: %s' % args.annotator)

    # elif args.annotator == Annotators.REFSEQ:
    #     annotator = RefSeq(
    #         gtf_file=args.refseq_gtf_file,
    #         version=args.refseq_version,
    #         species=args.refseq_species,
    #         assembly_report_txt_file=args.refseq_assembly_report_txt_file
    #     )

    # Step 3. Annotate variants
    variants_list = annotate(
        variants_list=variants_list,
        annotator=annotator,
        num_threads=args.num_threads
    )

    # Step 4. Write to a TSV file
    df_variants = variants_list.to_dataframe()
    df_variants.to_csv(
        args.output_tsv_file,
        sep='\t',
        index=False
    )
