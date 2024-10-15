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
and run 'filter' command.
"""


import argparse
from collections import defaultdict
from ..constants import VariantFilterSampleTypes, VariantCallTags
from ..default import *
from ..genomic_ranges_list import GenomicRangesList
from ..logging import get_logger
from ..main import filter, filter_excluded_regions, filter_homopolymeric_variants
from ..utilities import str2bool
from ..variant import Variant
from ..variant_filter import VariantFilter
from ..variants_list import VariantsList


logger = get_logger(__name__)


def add_cli_filter_arg_parser(
        sub_parsers: argparse._SubParsersAction
) -> argparse._SubParsersAction:
    """
    Adds 'filter' parser.
    """
    parser = sub_parsers.add_parser('filter', help='Filter a variants list.')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--tsv-file", '-i',
        dest="tsv_file",
        type=str,
        required=True,
        help="Input variants list TSV file."
    )
    parser_required.add_argument(
        "--case-sample-id", '-c',
        dest="case_sample_id",
        type=str,
        action='append',
        required=True,
        help="Case sample ID. Specify this parameter multiple times if there are "
             "multiple case sample IDs in the supplied TSV file. "
             "(e.g. --case-sample-id case_01 --case-sample-id case_02)."
    )
    parser_required.add_argument(
        "--output-passed-tsv-file", '-o',
        dest="output_passed_tsv_file",
        type=str,
        required=True,
        help="Output (passed) TSV file."
    )
    parser_required.add_argument(
        "--output-rejected-tsv-file", '-r',
        dest="output_rejected_tsv_file",
        type=str,
        required=True,
        help="Output (rejected) TSV file."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        '--control-sample-id',
        dest='control_sample_id',
        type=str,
        action='append',
        required=False,
        help="Control sample ID. Specify this parameter multiple times if there are "
             "multiple control sample IDs in the supplied TSV file. "
             "(e.g. --control-sample-id control_01 --control-sample-id control_02)."
    )
    parser_optional.add_argument(
        '--filter',
        dest='filter',
        type=str,
        action='append',
        required=False,
        help='Variant filter conditions: '
             '"{case,control} {all,average,median,min,max,any} {attribute} {<,<=,>,>=,==,in} {value}". '
             'Example 1: "case all alternate_allele_read_count >= 3". '
             'Example 2: "case all chr_1 in ["chr1","chr2","chr3"]". '
             'Please refer to the VSTOL documentation on how the filter semantics work.'
    )
    parser_optional.add_argument(
        "--reference-genome-fasta-file",
        dest="reference_genome_fasta_file",
        type=str,
        required=False,
        help="Reference genome FASTA file. If you specify this file, "
             "variant calls in homopolymer regions will be filtered out."
    )
    parser_optional.add_argument(
        "--excluded-regions-tsv-file",
        dest="excluded_regions_tsv_file",
        type=str,
        required=False,
        help="TSV files of regions to exclude. "
             "Variants with any variant calls where breakpoints are near the regions in this file will be removed. "
             "Expected headers: 'chromosome', 'start', 'end'."
    )
    parser_optional.add_argument(
        "--homopolymer-length",
        dest="homopolymer_length",
        type=int,
        required=False,
        default=HOMOPOLYMER_LENGTH,
        help="Variant calls with a homopolymer of this length in either breakpoint will be filtered out."
             "(default: %i)."
             % HOMOPOLYMER_LENGTH
    )
    parser_optional.add_argument(
        "--num-threads",
        dest="num_threads",
        type=int,
        default=NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)." % NUM_THREADS
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

    parser.set_defaults(which='filter')
    return sub_parsers


def run_cli_filter_from_parsed_args(args: argparse.Namespace):
    """
    Run 'filter' command using parameters from parsed arguments.

    Parameters:
        args    :   argparse.Namespace object with the following variables:
                    tsv_file
                    case_sample_id
                    output_filtered_tsv_file
                    output_rejected_tsv_file
                    control_sample_id
                    filter
                    reference_genome_fasta_file
                    excluded_regions_tsv_file
                    excluded_variants_tsv_file
                    excluded_variants_padding
                    homopolymer_length
                    num_threads
                    gzip
    """
    # Step 1. Load variants list
    logger.info('Loading target variants list')
    variants_list = VariantsList.read_tsv_file(tsv_file=args.tsv_file)

    # Step 2. Load variant filters
    variant_filters = []
    if args.filter is not None:
        for curr_filter in args.filter:
            curr_filter = curr_filter.split(' ')
            if curr_filter[0] == VariantFilterSampleTypes.CASE:
                sample_ids = args.case_sample_id
            elif curr_filter[0] == VariantFilterSampleTypes.CONTROL:
                sample_ids = args.control_sample_id
            else:
                raise Exception('Unknown sample type for variant filter: %s.' % curr_filter[0])
            variant_filter = VariantFilter(
                quantifier=curr_filter[1],
                attribute=curr_filter[2],
                operator=curr_filter[3],
                value=curr_filter[4],
                sample_ids=sample_ids
            )
            variant_filters.append(variant_filter)
    logger.info('Loaded %i variant filters.' % len(variant_filters))

    # Step 3. Load excluded regions list
    if args.excluded_regions_tsv_file is not None:
        excluded_regions_list = GenomicRangesList.read_tsv_file(tsv_file=args.excluded_regions_tsv_file)
        logger.info('Loaded %i excluded regions.' % excluded_regions_list.num_genomic_regions)
    else:
        excluded_regions_list = None

    # Step 4. Apply variant filters
    if args.filter is not None:
        _, variants_list_filters_rejected = filter(
            variants_list=variants_list,
            variant_filters=variant_filters,
            num_threads=args.num_threads
        )
    else:
        variants_list_filters_rejected = VariantsList()

    # Step 5. Filter out variant calls in excluded regions
    if excluded_regions_list is not None:
        _, variants_list_excluded_regions_rejected = filter_excluded_regions(
            variants_list=variants_list,
            excluded_regions_list=excluded_regions_list,
            num_threads=args.num_threads
        )
    else:
        variants_list_excluded_regions_rejected = VariantsList()

    # Step 6. Filter out homopolymeric variant calls
    if args.reference_genome_fasta_file is not None:
        _, variants_list_homopolymer_rejected = filter_homopolymeric_variants(
            variants_list=variants_list,
            reference_genome_fasta_file=args.reference_genome_fasta_file,
            homopolymer_length=args.homopolymer_length,
            num_threads=args.num_threads
        )
    else:
        variants_list_homopolymer_rejected = VariantsList()

    # Step 7. Consolidate the tags for rejected variant calls
    rejected_variant_call_ids_dict = defaultdict(list)
    for variant in variants_list_filters_rejected.variants:
        for variant_call in variant.variant_calls:
            rejected_variant_call_ids_dict[variant_call.id].append(VariantCallTags.FAILED_FILTER)
    for variant in variants_list_excluded_regions_rejected.variants:
        for variant_call in variant.variant_calls:
            rejected_variant_call_ids_dict[variant_call.id].append(VariantCallTags.NEARBY_EXCLUDED_REGION)
    for variant in variants_list_homopolymer_rejected.variants:
        for variant_call in variant.variant_calls:
            rejected_variant_call_ids_dict[variant_call.id].append(VariantCallTags.HOMOPOLYMER_REGION)

    # Step 8. Consolidate the passed and rejected variants lists
    variants_list_passed = VariantsList()
    variants_list_rejected = VariantsList()
    for variant in variants_list.variants:
        variant_calls_passed = []
        variant_calls_rejected = []
        for variant_call in variant.variant_calls:
            if variant_call.id in rejected_variant_call_ids_dict.keys():
                for tag in rejected_variant_call_ids_dict[variant_call.id]:
                    variant_call.tags.add(tag)
                variant_calls_rejected.append(variant_call)
            else:
                variant_call.tags.add(VariantCallTags.PASSED)
                variant_calls_passed.append(variant_call)
        variant_passed = Variant(id=variant.id)
        variant_rejected = Variant(id=variant.id)
        for variant_call in variant_calls_passed:
            variant_passed.add_variant_call(variant_call=variant_call)
        for variant_call in variant_calls_rejected:
            variant_rejected.add_variant_call(variant_call=variant_call)
        if variant_passed.num_variant_calls > 0:
            variants_list_passed.add_variant(variant=variant_passed)
        if variant_rejected.num_variant_calls > 0:
            variants_list_rejected.add_variant(variant=variant_rejected)

    # Step 9. Write to TSV files
    logger.info('Started writing passed and rejected variants lists')
    df_variants_passed = variants_list_passed.to_dataframe()
    df_variants_passed.sort_values(['variant_id'], inplace=True)
    df_variants_rejected = variants_list_rejected.to_dataframe()
    df_variants_rejected.sort_values(['variant_id'], inplace=True)
    if args.gzip:
        if args.output_passed_tsv_file.endswith(".gz") == False:
            args.output_passed_tsv_file = args.output_passed_tsv_file + '.gz'
        if args.output_rejected_tsv_file.endswith(".gz") == False:
            args.output_rejected_tsv_file = args.output_rejected_tsv_file + '.gz'
        df_variants_passed.to_csv(
            args.output_passed_tsv_file,
            sep='\t',
            index=False,
            compression='gzip'
        )
        df_variants_rejected.to_csv(
            args.output_rejected_tsv_file,
            sep='\t',
            index=False,
            compression='gzip'
        )
    else:
        df_variants_passed.to_csv(
            args.output_passed_tsv_file,
            sep='\t',
            index=False
        )
        df_variants_rejected.to_csv(
            args.output_rejected_tsv_file,
            sep='\t',
            index=False
        )
    logger.info('Finished writing passed and rejected variants lists')
