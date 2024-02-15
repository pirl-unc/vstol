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
The purpose of this python3 script is to implement the RefSeq dataclass.
"""


import gzip
import pandas as pd
import multiprocessing as mp
from collections import defaultdict
from dataclasses import dataclass
from functools import partial
from typing import Dict, List, Tuple
from .annotator import Annotator
from .constants import *
from .logging import get_logger
from .variant import Variant
from .variant_call import VariantCall
from .variant_call_annotation import VariantCallAnnotation
from .variants_list import VariantsList


logger = get_logger(__name__)


@dataclass
class RefSeq(Annotator):
    gtf_file: str
    assembly_report_txt_file: str
    version: str                            # RefSeq GTF file version (e.g. 'v110')
    species: str                            # RefSeq GTF file genome species (e.g. 'human')
    assembly_dict: Dict[str,str] = None
    df_genes: pd.DataFrame = None
    df_transcripts: pd.DataFrame = None
    df_exons: pd.DataFrame = None
    df_cds: pd.DataFrame = None
    df_start_codons: pd.DataFrame = None
    df_stop_codons: pd.DataFrame = None

    def __post_init__(self):
        self.__read_assembly_report_txt_file()
        self.__read_gtf_file_genes()              # add genes
        self.__read_gtf_file_transcripts()    # add transcripts
        self.__read_gtf_file_exons()              # add exons
        self.__read_gtf_file_cds()            # update UTR start and end positions
        self.__read_gtf_file_start_codons()
        self.__read_gtf_file_stop_codons()

    def __read_assembly_report_txt_file(self):
        """
        Read RefSeq assembly report TXT file and populate self.assembly_dict.
        """
        self.assembly_dict = {}
        with open(self.assembly_report_txt_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line[0:1] == '#':
                    continue
                elements = line.split('\t')
                genbank_accn = elements[6]
                ucsc_chromosome = elements[9]
                self.assembly_dict[genbank_accn] = ucsc_chromosome

    def __read_gtf_file_genes(self):
        """
        Read RefSeq GTF file and write rows to self.df_genes.
        """
        data = {
            'gene_id': [],
            'source': [],
            'name': [],
            'chromosome': [],
            'start': [],
            'end': [],
            'strand': [],
            'type': []
        }
        with open(self.gtf_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line[0:1] == '#':
                    continue
                elements = line.split('\t')
                if elements[2] == 'gene':
                    curr_gene_chrom = self.assembly_dict[str(elements[0])]
                    curr_gene_source = str(elements[1])
                    curr_gene_start = int(elements[3])
                    curr_gene_end = int(elements[4])
                    curr_gene_strand = str(elements[6])
                    curr_metadata = str(elements[8]).split(';')
                    curr_metadata_dict = {}
                    for curr_metadata_elements in curr_metadata:
                        if curr_metadata_elements == '':
                            continue
                        if curr_metadata_elements[0] == ' ':
                            curr_metadata_elements = curr_metadata_elements[1:]
                        curr_metadata_elements_ = curr_metadata_elements.split(' ')
                        curr_metadata_dict[curr_metadata_elements_[0]] = curr_metadata_elements_[1].replace('"', '')
                    curr_gene_id = str(curr_metadata_dict['gene'])
                    curr_gene_name = str(curr_metadata_dict['gene'])
                    curr_gene_type = str(curr_metadata_dict['gene_biotype'])
                    data['gene_id'].append(curr_gene_id)
                    data['source'].append(curr_gene_source)
                    data['name'].append(curr_gene_name)
                    data['chromosome'].append(curr_gene_chrom)
                    data['start'].append(curr_gene_start)
                    data['end'].append(curr_gene_end)
                    data['strand'].append(curr_gene_strand)
                    data['type'].append(curr_gene_type)
        self.df_genes = pd.DataFrame(data)
        logger.info('Loaded %i genes in total.' % len(self.df_genes))

    def __read_gtf_file_transcripts(self):
        """
        Read RefSeq GTF file and write rows to self.df_transcripts.
        """
        data = {
            'gene_id': [],
            'transcript_id': [],
            'source': [],
            'chromosome': [],
            'start': [],
            'end': [],
            'type': [],
            'strand': []
        }
        with open(self.gtf_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line[0:1] == '#':
                    continue
                elements = line.split('\t')
                if elements[2] == 'transcript':
                    curr_transcript_chrom = self.assembly_dict[str(elements[0])]
                    curr_transcript_source = str(elements[1])
                    curr_transcript_start = int(elements[3])
                    curr_transcript_end = int(elements[4])
                    curr_transcript_strand = str(elements[6])
                    curr_metadata = str(elements[8]).split(';')
                    curr_metadata_dict = defaultdict(list)
                    for curr_metadata_elements in curr_metadata:
                        if curr_metadata_elements == '':
                            continue
                        if curr_metadata_elements[0] == ' ':
                            curr_metadata_elements = curr_metadata_elements[1:]
                        curr_metadata_elements_ = curr_metadata_elements.split(' ')
                        curr_metadata_dict[curr_metadata_elements_[0]].append(curr_metadata_elements_[1].replace('"', ''))
                    curr_gene_id = str(curr_metadata_dict['gene_id'][0])
                    curr_transcript_id = str(curr_metadata_dict['transcript_id'][0])
                    curr_transcript_type = str(curr_metadata_dict['transcript_biotype'][0])
                    data['gene_id'].append(curr_gene_id)
                    data['transcript_id'].append(curr_transcript_id)
                    data['source'].append(curr_transcript_source)
                    data['chromosome'].append(curr_transcript_chrom)
                    data['start'].append(curr_transcript_start)
                    data['end'].append(curr_transcript_end)
                    data['type'].append(curr_transcript_type)
                    data['strand'].append(curr_transcript_strand)
            self.df_transcripts = pd.DataFrame(data)
        logger.info('Loaded %i transcripts in total.' % len(self.df_transcripts))

    def __read_gtf_file_exons(self):
        """
        Read RefSeq GTF file and write rows to self.df_exons.
        """
        data = {
            'gene_id': [],
            'transcript_id': [],
            'exon_id': [],
            'source': [],
            'chromosome': [],
            'start': [],
            'end': [],
            'number': [],
            'strand': []
        }
        with open(self.gtf_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line[0:1] == '#':
                    continue
                elements = line.split('\t')
                if elements[2] == 'exon':
                    curr_exon_chrom = self.assembly_dict[str(elements[0])]
                    curr_exon_source = str(elements[1])
                    curr_exon_start = int(elements[3])
                    curr_exon_end = int(elements[4])
                    curr_exon_strand = str(elements[6])
                    curr_metadata = str(elements[8]).split(';')
                    curr_metadata_dict = defaultdict(list)
                    for curr_metadata_elements in curr_metadata:
                        if curr_metadata_elements == '':
                            continue
                        if curr_metadata_elements[0] == ' ':
                            curr_metadata_elements = curr_metadata_elements[1:]
                        curr_metadata_elements_ = curr_metadata_elements.split(' ')
                        curr_metadata_dict[curr_metadata_elements_[0]].append(curr_metadata_elements_[1].replace('"', ''))
                    curr_gene_id = str(curr_metadata_dict['gene_id'][0])
                    curr_transcript_id = str(curr_metadata_dict['transcript_id'][0])
                    curr_exon_id = '%s-%s-%s' % (curr_gene_id, curr_transcript_id, curr_metadata_dict['exon_number'][0])
                    curr_exon_number = int(curr_metadata_dict['exon_number'][0])
                    data['gene_id'].append(curr_gene_id)
                    data['transcript_id'].append(curr_transcript_id)
                    data['exon_id'].append(curr_exon_id)
                    data['source'].append(curr_exon_source)
                    data['chromosome'].append(curr_exon_chrom)
                    data['start'].append(curr_exon_start)
                    data['end'].append(curr_exon_end)
                    data['number'].append(curr_exon_number)
                    data['strand'].append(curr_exon_strand)
        self.df_exons = pd.DataFrame(data)
        logger.info('Loaded %i exons in total.' % len(self.df_exons))

    def __read_gtf_file_cds(self):
        """
        Read RefSeq GTF file and write rows to self.df_cds.
        """
        data = {
            'gene_id': [],
            'transcript_id': [],
            'exon_id': [],
            'cds_id': [],
            'source': [],
            'chromosome': [],
            'start': [],
            'end': [],
            'exon_number': [],
            'strand': []
        }
        with open(self.gtf_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line[0:1] == '#':
                    continue
                elements = line.split('\t')
                if elements[2] == 'CDS':
                    curr_cds_chrom = self.assembly_dict[str(elements[0])]
                    curr_cds_source = str(elements[1])
                    curr_cds_start = int(elements[3])
                    curr_cds_end = int(elements[4])
                    curr_cds_strand = str(elements[6])
                    curr_metadata = str(elements[8]).split(';')
                    curr_metadata_dict = defaultdict(list)
                    for curr_metadata_elements in curr_metadata:
                        if curr_metadata_elements == '':
                            continue
                        if curr_metadata_elements[0] == ' ':
                            curr_metadata_elements = curr_metadata_elements[1:]
                        curr_metadata_elements_ = curr_metadata_elements.split(' ')
                        curr_metadata_dict[curr_metadata_elements_[0]].append(curr_metadata_elements_[1].replace('"', ''))
                    curr_gene_id = str(curr_metadata_dict['gene_id'][0])
                    curr_transcript_id = str(curr_metadata_dict['transcript_id'][0])
                    curr_exon_id = int(curr_metadata_dict['exon_number'][0])
                    curr_exon_number = curr_exon_id
                    curr_cds_id = '%s-%s-%i' % (curr_gene_id, curr_transcript_id, curr_exon_number)
                    data['gene_id'].append(curr_gene_id)
                    data['transcript_id'].append(curr_transcript_id)
                    data['exon_id'].append(curr_exon_id)
                    data['cds_id'].append(curr_cds_id)
                    data['source'].append(curr_cds_source)
                    data['chromosome'].append(curr_cds_chrom)
                    data['start'].append(curr_cds_start)
                    data['end'].append(curr_cds_end)
                    data['exon_number'].append(curr_exon_number)
                    data['strand'].append(curr_cds_strand)
        self.df_cds = pd.DataFrame(data)
        logger.info('Loaded %i CDS in total.' % len(self.df_cds))

    def __read_gtf_file_start_codons(self):
        """
        Read RefSeq GTF file and write rows to self.df_start_codons.
        """
        data = {
            'gene_id': [],
            'transcript_id': [],
            'exon_id': [],
            'start_codon_id': [],
            'source': [],
            'chromosome': [],
            'start': [],
            'end': [],
            'exon_number': [],
            'strand': []
        }
        with open(self.gtf_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line[0:1] == '#':
                    continue
                elements = line.split('\t')
                if elements[2] == 'start_codon':
                    curr_start_codon_chrom = self.assembly_dict[str(elements[0])]
                    curr_start_codon_source = str(elements[1])
                    curr_start_codon_start = int(elements[3])
                    curr_start_codon_end = int(elements[4])
                    curr_start_codon_strand = str(elements[6])
                    curr_metadata = str(elements[8]).split(';')
                    curr_metadata_dict = defaultdict(list)
                    for curr_metadata_elements in curr_metadata:
                        if curr_metadata_elements == '':
                            continue
                        if curr_metadata_elements[0] == ' ':
                            curr_metadata_elements = curr_metadata_elements[1:]
                        curr_metadata_elements_ = curr_metadata_elements.split(' ')
                        curr_metadata_dict[curr_metadata_elements_[0]].append(curr_metadata_elements_[1].replace('"', ''))
                    curr_gene_id = str(curr_metadata_dict['gene_id'][0])
                    curr_transcript_id = str(curr_metadata_dict['transcript_id'][0])
                    curr_exon_id = int(curr_metadata_dict['exon_number'][0])
                    curr_exon_number = curr_exon_id
                    curr_start_codon_id = '%s-%s-%i' % (curr_gene_id, curr_transcript_id, curr_exon_number)
                    data['gene_id'].append(curr_gene_id)
                    data['transcript_id'].append(curr_transcript_id)
                    data['exon_id'].append(curr_exon_id)
                    data['start_codon_id'].append(curr_start_codon_id)
                    data['source'].append(curr_start_codon_source)
                    data['chromosome'].append(curr_start_codon_chrom)
                    data['start'].append(curr_start_codon_start)
                    data['end'].append(curr_start_codon_end)
                    data['exon_number'].append(curr_exon_number)
                    data['strand'].append(curr_start_codon_strand)
        self.df_start_codons = pd.DataFrame(data)
        logger.info('Loaded %i start codons in total.' % len(self.df_start_codons))

    def __read_gtf_file_stop_codons(self):
        """
        Read RefSeq GTF file and write rows to self.df_start_codons.
        """
        data = {
            'gene_id': [],
            'transcript_id': [],
            'exon_id': [],
            'stop_codon_id': [],
            'source': [],
            'chromosome': [],
            'start': [],
            'end': [],
            'exon_number': [],
            'strand': []
        }
        with open(self.gtf_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line[0:1] == '#':
                    continue
                elements = line.split('\t')
                if elements[2] == 'stop_codon':
                    curr_stop_codon_chrom = self.assembly_dict[str(elements[0])]
                    curr_stop_codon_source = str(elements[1])
                    curr_stop_codon_start = int(elements[3])
                    curr_stop_codon_end = int(elements[4])
                    curr_stop_codon_strand = str(elements[6])
                    curr_metadata = str(elements[8]).split(';')
                    curr_metadata_dict = defaultdict(list)
                    for curr_metadata_elements in curr_metadata:
                        if curr_metadata_elements == '':
                            continue
                        if curr_metadata_elements[0] == ' ':
                            curr_metadata_elements = curr_metadata_elements[1:]
                        curr_metadata_elements_ = curr_metadata_elements.split(' ')
                        curr_metadata_dict[curr_metadata_elements_[0]].append(curr_metadata_elements_[1].replace('"', ''))
                    curr_gene_id = str(curr_metadata_dict['gene_id'][0])
                    curr_transcript_id = str(curr_metadata_dict['transcript_id'][0])
                    curr_exon_id = int(curr_metadata_dict['exon_number'][0])
                    curr_exon_number = curr_exon_id
                    curr_stop_codon_id = '%s-%s-%i' % (curr_gene_id, curr_transcript_id, curr_exon_number)
                    data['gene_id'].append(curr_gene_id)
                    data['transcript_id'].append(curr_transcript_id)
                    data['exon_id'].append(curr_exon_id)
                    data['stop_codon_id'].append(curr_stop_codon_id)
                    data['source'].append(curr_stop_codon_source)
                    data['chromosome'].append(curr_stop_codon_chrom)
                    data['start'].append(curr_stop_codon_start)
                    data['end'].append(curr_stop_codon_end)
                    data['exon_number'].append(curr_exon_number)
                    data['strand'].append(curr_stop_codon_strand)
        self.df_stop_codons = pd.DataFrame(data)
        logger.info('Loaded %i stop codons in total.' % len(self.df_stop_codons))

    def annotate_position(
            self,
            chromosome: str,
            position: int
    ) -> List[VariantCallAnnotation]:
        """
        Annotate a position and return a list of VariantCallAnnotation objects.

        Parameters:
            chromosome              :   Chromosome.
            position                :   Position.

        Returns:
            List[VariantCallAnnotation]
        """
        variant_call_annotations = []
        df_genes_matched = self.df_genes[
            (self.df_genes['chromosome'] == chromosome) &
            (self.df_genes['start'] <= position) &
            (self.df_genes['end'] >= position)
        ]
        if len(df_genes_matched) == 0:
            variant_call_annotation = VariantCallAnnotation(
                annotator=Annotators.REFSEQ,
                annotator_version=self.version,
                region=GenomicRegionTypes.INTERGENIC,
                species=self.species
            )
            variant_call_annotations.append(variant_call_annotation)
        else:
            for _, row in df_genes_matched.iterrows():
                df_cds_matched = self.df_cds[
                    (self.df_cds['gene_id'] == row['gene_id']) &
                    (self.df_cds['start'] <= position) &
                    (self.df_cds['end'] >= position)
                ]
                df_exons_matched = self.df_exons[
                    (self.df_exons['gene_id'] == row['gene_id']) &
                    (self.df_exons['start'] <= position) &
                    (self.df_exons['end'] >= position)
                ]
                df_start_codon_matched = self.df_start_codons[
                    (self.df_start_codons['gene_id'] == row['gene_id']) &
                    (self.df_start_codons['start'] <= position) &
                    (self.df_start_codons['end'] >= position)
                ]
                df_stop_codon_matched = self.df_stop_codons[
                    (self.df_stop_codons['gene_id'] == row['gene_id']) &
                    (self.df_stop_codons['start'] <= position) &
                    (self.df_stop_codons['end'] >= position)
                ]
                if len(df_cds_matched) > 0:
                    if len(df_start_codon_matched) > 0:
                        region = GenomicRegionTypes.START_CODON
                    elif len(df_stop_codon_matched) > 0:
                        region = GenomicRegionTypes.STOP_CODON
                    else:
                        region = GenomicRegionTypes.EXONIC
                else:
                    if len(df_exons_matched) > 0:
                        region = GenomicRegionTypes.UNTRANSLATED_REGION
                    else:
                        region = GenomicRegionTypes.INTRONIC
                variant_call_annotation = VariantCallAnnotation(
                    annotator=Annotators.REFSEQ,
                    annotator_version=self.version,
                    gene_id=row['gene_id'],
                    gene_id_stable=row['gene_id'],
                    gene_name=row['name'],
                    gene_strand=row['strand'],
                    gene_type=row['type'],
                    region=region,
                    species=self.species
                )
                variant_call_annotations.append(variant_call_annotation)
        return variant_call_annotations

    def annotate_variant_call(
            self,
            variant_call: VariantCall
    ) -> Tuple[List[VariantCallAnnotation], List[VariantCallAnnotation]]:
        """
        Annotate a VariantCall object and return two lists of VariantCallAnnotation objects.

        Parameters:
            variant_call        :   VariantCall object.

        Returns:
            Tuple[pos_1_annotations,pos_2_annotations]
        """
        pos_1_annotations = self.annotate_position(
            chromosome=variant_call.chromosome_1,
            position=variant_call.position_1
        )
        pos_2_annotations = self.annotate_position(
            chromosome=variant_call.chromosome_2,
            position=variant_call.position_2
        )
        return pos_1_annotations, pos_2_annotations

    def annotate_variant(self, variant: Variant) -> Variant:
        for i in range(0, variant.num_variant_calls):
            position_1_annotations, position_2_annotations = self.annotate_variant_call(
                variant.variant_calls[i]
            )
            for annotation in position_1_annotations:
                variant.variant_calls[i].position_1_annotations.append(annotation)
            for annotation in position_2_annotations:
                variant.variant_calls[i].position_2_annotations.append(annotation)
        return variant

    def annotate(self,
                 variants_list,
                 num_processes: int = 1) -> VariantsList:
        """
        Annotate a VariantsList object.

        Parameters:
            variants_list   :   VariantsList.
            num_processes   :   Number of processes.

        Returns:
            VariantsList
        """
        pool = mp.Pool(processes=num_processes)
        func = partial(self.annotate_variant)
        variants = pool.map(func, variants_list.variants)
        pool.close()
        variants_list_annotated = VariantsList()
        for variant in variants:
            variants_list_annotated.add_variant(variant=variant)
        return variants_list_annotated
