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
The purpose of this python3 script is to implement a parser for GRIDSS2 VCF files.
"""


from .common import *
from ..constants import *
from ..logging import get_logger
from ..utilities import retrieve_from_dict


logger = get_logger(__name__)


def parse_gridss2_vcf(
        df_vcf: pd.DataFrame,
        source_id: str,
        platform: str
) -> pd.DataFrame:
    """
    Parse a GRIDSS2 VCF DataFrame and return a Pandas DataFrame.
    """
    # Step 1. Get a template dictionary
    data = get_template_dict()

    # Step 2. Get all sample IDs
    sample_ids = df_vcf.columns.values.tolist()[9:]

    def info_to_dict(info_list):
        d = {}
        for curr_info in info_list:
            if '=' in curr_info:
                k, v = curr_info.split('=', 1)
                d[k] = str(v)
            else:
                d[curr_info] = True
        return d

    # Step 3. Record variant calls
    per_sample_events = {sid: {} for sid in sample_ids}
    for row in df_vcf.to_dict('records'):
        chrom = retrieve_from_dict(dct=row, key='CHROM', default_value='', type=str)
        pos   = retrieve_from_dict(dct=row, key='POS', default_value='', type=int)
        ref   = retrieve_from_dict(dct=row, key='REF', default_value='', type=str)
        alt   = retrieve_from_dict(dct=row, key='ALT', default_value='', type=str)
        qual  = retrieve_from_dict(dct=row, key='QUAL', default_value=-1.0, type=float)
        filter  = retrieve_from_dict(dct=row, key='FILTER', default_value='', type=str)

        info_list = str(row['INFO']).split(';')
        fmt_keys  = str(row['FORMAT']).split(':')
        info_dict = info_to_dict(info_list)

        sv_type = info_dict.get('SVTYPE', None)
        if sv_type != 'BND':
            continue

        event_id = info_dict.get('EVENT', row.get('ID', ''))
        if event_id == '':
            # Fallback: use CHROM:POS:ALT
            event_id = f"{chrom}:{pos}:{alt}"

        for sample_id in sample_ids:
            fmt_vals = str(row[sample_id]).split(':')
            format_dict = {k: (fmt_vals[i] if i < len(fmt_vals) else '.') for i, k in enumerate(fmt_keys)}

            af = float(format_dict['AF']) if format_dict.get('AF', '.') != '.' else -1.0
            num_ref_reads = int(format_dict['REF']) if format_dict.get('REF', '.') != '.' else -1
            num_alt_reads = int(format_dict['VF']) if format_dict.get('VF', '.') != '.' else -1
            num_total_reads = (num_ref_reads + num_alt_reads) if (num_ref_reads != -1 and num_alt_reads != -1) else -1

            alt_bp = parse_alt_breakpoints(alternate_allele=alt, chromosome_1=chrom, position_1=pos)

            per_sample_events[sample_id].setdefault(event_id, []).append({
                "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                "qual": qual, "filter": filter,
                "info_dict": info_dict,
                "format_dict": format_dict,
                "af": af,
                "num_ref_reads": num_ref_reads,
                "num_alt_reads": num_alt_reads,
                "num_total_reads": num_total_reads,
                "alt_bp": alt_bp
            })

    variant_id = 1
    def append_variant(
            sample_id,
            qual,
            filter,
            info_dict,
            format_dict,
            chromosome_1,
            position_1,
            strand_1,
            operation_1,
            chromosome_2,
            position_2,
            strand_2,
            operation_2,
            sequence,
            size,
            variant_type,
            num_ref_reads,
            num_alt_reads,
            num_total_reads,
            af
    ):
        data['variant_id'].append(variant_id)
        data['source_id'].append(source_id)
        data['sample_id'].append(sample_id)
        data['platform'].append(platform)
        data['method'].append(VariantCallingMethod.GRIDSS2)
        data['quality'].append(qual)
        data['filter'].append(filter)
        data['chromosome_1'].append(chromosome_1)
        data['position_1'].append(position_1)
        data['strand_1'].append(strand_1)
        data['operation_1'].append(operation_1)
        data['chromosome_2'].append(chromosome_2)
        data['position_2'].append(position_2)
        data['strand_2'].append(strand_2)
        data['operation_2'].append(operation_2)
        data['sequence'].append(sequence)
        data['variant_size'].append(size)
        data['variant_type'].append(variant_type)
        data['info'].append(';'.join('%s=%s' % (k, v) for k, v in info_dict.items()))
        data['format'].append(';'.join('%s=%s' % (k, v) for k, v in format_dict.items()))
        data['position_1_read_count_reference_allele'].append(num_ref_reads)
        data['position_2_read_count_reference_allele'].append(num_ref_reads)
        data['position_1_read_count_alternate_allele'].append(num_alt_reads)
        data['position_2_read_count_alternate_allele'].append(num_alt_reads)
        data['position_1_read_count_total'].append(num_total_reads)
        data['position_2_read_count_total'].append(num_total_reads)
        data['position_1_alternate_allele_fraction'].append(af)
        data['position_2_alternate_allele_fraction'].append(af)
        data['read_ids_alternate_allele'].append('')

    for sample_id, event_map in per_sample_events.items():
        for event_id, recs in event_map.items():
            # Common: 2 records per event; sometimes 1 (single breakend) or >2 (complex)
            if len(recs) == 1:
                r = recs[0]
                bp = r["alt_bp"]
                append_variant(
                    sample_id=sample_id,
                    qual=r["qual"],
                    filter=r["filter"],
                    info_dict=r["info_dict"],
                    format_dict=r["format_dict"],
                    chromosome_1=bp["chromosome_1"],
                    position_1=bp["position_1"],
                    strand_1=bp["strand_1"],
                    operation_1=bp["operation_1"],
                    chromosome_2=bp["chromosome_2"],
                    position_2=bp["position_2"],
                    strand_2=bp["strand_2"],
                    operation_2=bp["operation_2"],
                    sequence=bp["sequence"],
                    size=-1,
                    variant_type=VariantType.BREAKPOINT,   # single-ended breakpoint
                    num_ref_reads=r["num_ref_reads"],
                    num_alt_reads=r["num_alt_reads"],
                    num_total_reads=r["num_total_reads"],
                    af=r["af"]
                )
                variant_id += 1
                continue

            # Choose two representative ends for simple INS/DEL detection
            # (for complex events, fall back to emitting each breakpoint)
            r1, r2 = recs[0], recs[1]
            bp1, bp2 = r1["alt_bp"], r2["alt_bp"]

            # If either mate is unknown (single-breakend in group), treat as breakpoints
            if bp1["position_2"] == -1 or bp2["position_2"] == -1 or bp1["chromosome_2"] == "." or bp2["chromosome_2"] == ".":
                for r in recs:
                    bp = r["alt_bp"]
                    append_variant(
                        sample_id=sample_id,
                        qual=r["qual"],
                        filter=r["filter"],
                        info_dict=r["info_dict"],
                        format_dict=r["format_dict"],
                        chromosome_1=bp["chromosome_1"],
                        position_1=bp["position_1"],
                        strand_1=bp["strand_1"],
                        operation_1=bp["operation_1"],
                        chromosome_2=bp["chromosome_2"],
                        position_2=bp["position_2"],
                        strand_2=bp["strand_2"],
                        operation_2=bp["operation_2"],
                        sequence=bp["sequence"],
                        size=-1,
                        variant_type=VariantType.BREAKPOINT,
                        num_ref_reads=r["num_ref_reads"], num_alt_reads=r["num_alt_reads"],
                        num_total_reads=r["num_total_reads"], af=r["af"]
                    )
                    variant_id += 1
                continue

            # Simple same-chromosome INS/DEL detection using the resolved mate from bp1
            chrA = bp1["chromosome_1"]
            chrB = bp1["chromosome_2"]
            posA = bp1["position_1"]
            posB = bp1["position_2"]

            if chrA == chrB:
                lo, hi = (posA, posB) if posA <= posB else (posB, posA)
                dist = hi - lo

                # Insertion: adjacent breakpoints (pos2 == pos1 + 1)
                if dist == 1:
                    ins_seq = bp1.get("sequence", "") or bp2.get("sequence", "")
                    append_variant(
                        sample_id=sample_id,
                        qual=r1["qual"],
                        filter=r1["filter"],
                        info_dict=r1["info_dict"],
                        format_dict=r1["format_dict"],
                        chromosome_1=chrA,
                        position_1=lo,
                        strand_1=Strand.FORWARD,
                        operation_1=OperationType.DOWNSTREAM,
                        chromosome_2=chrA,
                        position_2=hi,
                        strand_2=Strand.FORWARD,
                        operation_2=OperationType.UPSTREAM,
                        sequence=ins_seq,
                        size=len(ins_seq),
                        variant_type=VariantType.INSERTION if hasattr(VariantType, "INSERTION") else VariantType.BREAKPOINT,
                        num_ref_reads=r1["num_ref_reads"],
                        num_alt_reads=r1["num_alt_reads"],
                        num_total_reads=r1["num_total_reads"],
                        af=r1["af"]
                    )
                    variant_id += 1
                    continue

                # Deletion: gap between breakpoints (pos2 > pos1 + 1)
                if dist > 1:
                    del_size = dist - 1
                    append_variant(
                        sample_id=sample_id,
                        qual=r1["qual"],
                        filter=r1["filter"],
                        info_dict=r1["info_dict"],
                        format_dict=r1["format_dict"],
                        chromosome_1=chrA,
                        position_1=lo,
                        strand_1=Strand.FORWARD,
                        operation_1=OperationType.DOWNSTREAM,
                        chromosome_2=chrA,
                        position_2=hi,
                        strand_2=Strand.FORWARD,
                        operation_2=OperationType.UPSTREAM,
                        sequence="",
                        size=del_size,
                        variant_type=VariantType.DELETION if hasattr(VariantType, "DELETION") else VariantType.BREAKPOINT,
                        num_ref_reads=r1["num_ref_reads"],
                        num_alt_reads=r1["num_alt_reads"],
                        num_total_reads=r1["num_total_reads"],
                        af=r1["af"]
                    )
                    variant_id += 1
                    continue

            # Otherwise, emit as breakpoint (TRA/INV/complex)
            for r in recs:
                bp = r["alt_bp"]
                append_variant(
                    sample_id=sample_id,
                    qual=r["qual"],
                    filter=r["filter"],
                    info_dict=r["info_dict"],
                    format_dict=r["format_dict"],
                    chromosome_1=bp["chromosome_1"],
                    position_1=bp["position_1"],
                    strand_1=bp["strand_1"],
                    operation_1=bp["operation_1"],
                    chromosome_2=bp["chromosome_2"],
                    position_2=bp["position_2"],
                    strand_2=bp["strand_2"],
                    operation_2=bp["operation_2"],
                    sequence=bp["sequence"],
                    size=abs(bp["position_1"] - bp["position_2"]) if bp["position_2"] != -1 else -1,
                    variant_type=VariantType.BREAKPOINT,
                    num_ref_reads=r["num_ref_reads"],
                    num_alt_reads=r["num_alt_reads"],
                    num_total_reads=r["num_total_reads"],
                    af=r["af"]
                )
                variant_id += 1

    df_variants = pd.DataFrame(data)
    logger.info('Converted %i GRIDSS2 variants.' % len(df_variants['variant_id'].unique()))
    return df_variants
