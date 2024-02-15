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
The purpose of this python3 script is to implement the Gencode dataclass.
"""


import pandas as pd
import multiprocessing as mp
from collections import defaultdict
from dataclasses import dataclass, field
from functools import partial
from typing import List, Tuple, Set
from .annotator import Annotator
from .constants import *
from .logging import get_logger
from .variant import Variant
from .variant_call import VariantCall
from .variant_call_annotation import VariantCallAnnotation
from .variants_list import VariantsList


logger = get_logger(__name__)


@dataclass
class Gencode(Annotator):
    gtf_file: str
    version: str    # GENCODE GTF file version (e.g. 'v41')
    species: str    # GENCODE GTF file genome species (e.g. 'human')
    levels: List[int] = field(default_factory=list)
    types: List[str] = field(default_factory=list)
    df_genes: pd.DataFrame = None
    df_transcripts: pd.DataFrame = None
    df_exons: pd.DataFrame = None
    df_start_codons: pd.DataFrame = None
    df_stop_codons: pd.DataFrame = None
    df_utrs: pd.DataFrame = None

    def __post_init__(self):
        self.__read_gtf_file_genes()  # add genes
        self.__read_gtf_file_transcripts()  # add transcripts
        self.__read_gtf_file_exons()  # add exons
        self.__read_gtf_file_start_codons()  # update start codon start and end positions
        self.__read_gtf_file_stop_codons()  # update stop codon start and end positions
        self.__read_gtf_file_utr()  # update UTR start and end positions

    @staticmethod
    def get_stable_ensembl_id(id: str) -> Tuple[str,str]:
        """
        Returns the stable Ensembl ID and version.

        Parameters:
            id          :   Ensembl ID (e.g. 'ENSG00001.1').

        Returns:
            Tuple[stable_id,version]
        """
        if (('ENSG' in id) or ('ENST' in id or 'ENSE' in id)) and '.' in id:
            stable_id = id.split('.')[0]
            version = id.split('.')[1]
            return stable_id, version
        else:
            return id, ''

    def __read_gtf_file_genes(self):
        """
        Read GENCODE GTF file and write rows to self.df_genes.
        """
        data = {
            'gene_id': [],
            'gene_id_stable': [],
            'source': [],
            'name': [],
            'chromosome': [],
            'start': [],
            'end': [],
            'strand': [],
            'type': [],
            'level': [],
            'version': []
        }
        with open(self.gtf_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line[0:2] == '##':
                    continue
                elements = line.split('\t')
                if elements[2] == 'gene':
                    curr_gene_chrom = str(elements[0])
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
                    curr_gene_id = str(curr_metadata_dict['gene_id'])
                    curr_gene_stable_id, curr_gene_version = Gencode.get_stable_ensembl_id(id=str(curr_metadata_dict['gene_id']))
                    curr_gene_name = str(curr_metadata_dict['gene_name'])
                    curr_gene_type = str(curr_metadata_dict['gene_type'])
                    curr_gene_level = int(curr_metadata_dict['level'])

                    if len(self.types) > 0 and curr_gene_type not in self.types:
                        continue
                    if len(self.levels) > 0 and curr_gene_level not in self.levels:
                        continue

                    data['gene_id'].append(curr_gene_id)
                    data['gene_id_stable'].append(curr_gene_stable_id)
                    data['source'].append(curr_gene_source)
                    data['name'].append(curr_gene_name)
                    data['chromosome'].append(curr_gene_chrom)
                    data['start'].append(curr_gene_start)
                    data['end'].append(curr_gene_end)
                    data['strand'].append(curr_gene_strand)
                    data['type'].append(curr_gene_type)
                    data['level'].append(curr_gene_level)
                    data['version'].append(curr_gene_version)
            self.df_genes = pd.DataFrame(data)
        self.__gene_ids = set()
        for gene_id in self.df_genes['gene_id'].values.tolist():
            self.__gene_ids.add(gene_id)
        logger.info('Loaded %i genes in total.' % len(self.df_genes))

    def __read_gtf_file_transcripts(self):
        """
        Read GENCODE GTF file and write rows to self.df_genes.
        """
        data = {
            'gene_id': [],
            'transcript_id': [],
            'transcript_id_stable': [],
            'source': [],
            'chromosome': [],
            'start': [],
            'end': [],
            'type': [],
            'strand': [],
            'version': [],
            'name': [],
            'level': [],
            'support_level': [],
            'tags': []
        }
        with open(self.gtf_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line[0:2] == '##':
                    continue
                elements = line.split('\t')
                if elements[2] == 'transcript':
                    curr_transcript_chrom = str(elements[0])
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
                    curr_transcript_stable_id, curr_transcript_version = Gencode.get_stable_ensembl_id(id=str(curr_metadata_dict['transcript_id'][0]))
                    curr_transcript_type = str(curr_metadata_dict['transcript_type'][0])
                    curr_transcript_name = str(curr_metadata_dict['transcript_name'][0])
                    curr_transcript_tags = [str(tag).replace('"', '') for tag in curr_metadata_dict['tag']]
                    try:
                        curr_transcript_level = int(curr_metadata_dict['level'][0])
                    except:
                        curr_transcript_level = ''
                    try:
                        curr_transcript_support_level = int(curr_metadata_dict['transcript_support_level'][0])
                    except:
                        curr_transcript_support_level = ''
                    if curr_gene_id not in self.__gene_ids:
                        continue
                    if len(self.types) > 0 and curr_transcript_type not in self.types:
                        continue
                    if len(self.levels) > 0 and curr_transcript_level not in self.levels:
                        continue
                    data['gene_id'].append(curr_gene_id)
                    data['transcript_id'].append(curr_transcript_id)
                    data['transcript_id_stable'].append(curr_transcript_stable_id)
                    data['source'].append(curr_transcript_source)
                    data['chromosome'].append(curr_transcript_chrom)
                    data['start'].append(curr_transcript_start)
                    data['end'].append(curr_transcript_end)
                    data['type'].append(curr_transcript_type)
                    data['strand'].append(curr_transcript_strand)
                    data['version'].append(curr_transcript_version)
                    data['name'].append(curr_transcript_name)
                    data['level'].append(curr_transcript_level)
                    data['support_level'].append(curr_transcript_support_level)
                    data['tags'].append(';'.join(curr_transcript_tags))
            self.df_transcripts = pd.DataFrame(data)
        self.__transcript_ids = set()
        for transcript_id in self.df_transcripts['transcript_id'].values.tolist():
            self.__transcript_ids.add(transcript_id)
        logger.info('Loaded %i transcripts in total.' % len(self.df_transcripts))

    def __read_gtf_file_exons(self):
        """
        Read GENCODE GTF file and write rows to self.df_exons.
        """
        data = {
            'gene_id': [],
            'transcript_id': [],
            'exon_id': [],
            'exon_id_stable': [],
            'source': [],
            'chromosome': [],
            'start': [],
            'end': [],
            'number': [],
            'strand': [],
            'version': [],
            'tags': []
        }
        with open(self.gtf_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line[0:2] == '##':
                    continue
                elements = line.split('\t')
                if elements[2] == 'exon':
                    curr_exon_chrom = str(elements[0])
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
                    curr_exon_id = str(curr_metadata_dict['exon_id'][0])
                    curr_exon_stable_id, curr_exon_version = Gencode.get_stable_ensembl_id(id=str(curr_metadata_dict['exon_id'][0]))
                    curr_exon_number = int(curr_metadata_dict['exon_number'][0])
                    curr_exon_tags = [str(tag).replace('"', '') for tag in curr_metadata_dict['tag']]
                    if curr_gene_id not in self.__gene_ids:
                        continue
                    if curr_transcript_id not in self.__transcript_ids:
                        continue
                    data['gene_id'].append(curr_gene_id)
                    data['transcript_id'].append(curr_transcript_id)
                    data['exon_id'].append(curr_exon_id)
                    data['exon_id_stable'].append(curr_exon_stable_id)
                    data['source'].append(curr_exon_source)
                    data['chromosome'].append(curr_exon_chrom)
                    data['start'].append(curr_exon_start)
                    data['end'].append(curr_exon_end)
                    data['number'].append(curr_exon_number)
                    data['strand'].append(curr_exon_strand)
                    data['version'].append(curr_exon_version)
                    data['tags'].append(';'.join(curr_exon_tags))
        self.df_exons = pd.DataFrame(data)
        logger.info('Loaded %i exons in total.' % len(self.df_exons))

    def __read_gtf_file_start_codons(self):
        """
        Read GENCODE GTF file and write rows to self.df_start_codons.
        """
        data = {
            'gene_id': [],
            'transcript_id': [],
            'start_codon_start': [],
            'start_codon_end': []
        }
        with open(self.gtf_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line[0:2] == '##':
                    continue
                elements = line.split('\t')
                if elements[2] == 'start_codon':
                    curr_start_codon_start = int(elements[3])
                    curr_start_codon_end = int(elements[4])
                    curr_metadata = str(elements[8]).split(';')
                    curr_metadata_dict = {}
                    for curr_metadata_elements in curr_metadata:
                        if curr_metadata_elements == '':
                            continue
                        if curr_metadata_elements[0] == ' ':
                            curr_metadata_elements = curr_metadata_elements[1:]
                        curr_metadata_elements_ = curr_metadata_elements.split(' ')
                        curr_metadata_dict[curr_metadata_elements_[0]] = curr_metadata_elements_[1].replace('"', '')
                    curr_gene_id = str(curr_metadata_dict['gene_id'])
                    curr_transcript_id = str(curr_metadata_dict['transcript_id'])
                    if curr_gene_id not in self.__gene_ids:
                        continue
                    if curr_transcript_id not in self.__transcript_ids:
                        continue
                    data['gene_id'].append(curr_gene_id)
                    data['transcript_id'].append(curr_transcript_id)
                    data['start_codon_start'].append(curr_start_codon_start)
                    data['start_codon_end'].append(curr_start_codon_end)
        self.df_start_codons = pd.DataFrame(data)
        logger.info('Loaded %i start codons in total.' % len(self.df_start_codons))

    def __read_gtf_file_stop_codons(self):
        """
        Read GENCODE GTF file and write rows to self.df_stop_codons.
        """
        data = {
            'gene_id': [],
            'transcript_id': [],
            'stop_codon_start': [],
            'stop_codon_end': []
        }
        with open(self.gtf_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line[0:2] == '##':
                    continue
                elements = line.split('\t')
                if elements[2] == 'stop_codon':
                    curr_stop_codon_start = int(elements[3])
                    curr_stop_codon_end = int(elements[4])
                    curr_metadata = str(elements[8]).split(';')
                    curr_metadata_dict = {}
                    for curr_metadata_elements in curr_metadata:
                        if curr_metadata_elements == '':
                            continue
                        if curr_metadata_elements[0] == ' ':
                            curr_metadata_elements = curr_metadata_elements[1:]
                        curr_metadata_elements_ = curr_metadata_elements.split(' ')
                        curr_metadata_dict[curr_metadata_elements_[0]] = curr_metadata_elements_[1].replace('"', '')
                    curr_gene_id = str(curr_metadata_dict['gene_id'])
                    curr_transcript_id = str(curr_metadata_dict['transcript_id'])
                    if curr_gene_id not in self.__gene_ids:
                        continue
                    if curr_transcript_id not in self.__transcript_ids:
                        continue
                    data['gene_id'].append(curr_gene_id)
                    data['transcript_id'].append(curr_transcript_id)
                    data['stop_codon_start'].append(curr_stop_codon_start)
                    data['stop_codon_end'].append(curr_stop_codon_end)
        self.df_stop_codons = pd.DataFrame(data)
        logger.info('Loaded %i stop codons in total.' % len(self.df_stop_codons))

    def __read_gtf_file_utr(self):
        """
        Read GENCODE GTF file and write rows to self.df_utrs.
        """
        data = {
            'gene_id': [],
            'transcript_id': [],
            'utr_start': [],
            'utr_end': [],
            'utr_type': []
        }
        with open(self.gtf_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line[0:2] == '##':
                    continue
                elements = line.split('\t')
                if elements[2] == 'UTR':
                    curr_utr_start = int(elements[3])
                    curr_utr_end = int(elements[4])
                    curr_metadata = str(elements[8]).split(';')
                    curr_metadata_dict = {}
                    for curr_metadata_elements in curr_metadata:
                        if curr_metadata_elements == '':
                            continue
                        if curr_metadata_elements[0] == ' ':
                            curr_metadata_elements = curr_metadata_elements[1:]
                        curr_metadata_elements_ = curr_metadata_elements.split(' ')
                        curr_metadata_dict[curr_metadata_elements_[0]] = curr_metadata_elements_[1].replace('"', '')
                    curr_gene_id = str(curr_metadata_dict['gene_id'])
                    curr_transcript_id = str(curr_metadata_dict['transcript_id'])
                    if curr_gene_id not in self.__gene_ids:
                        continue
                    if curr_transcript_id not in self.__transcript_ids:
                        continue
                    data['gene_id'].append(curr_gene_id)
                    data['transcript_id'].append(curr_transcript_id)
                    data['utr_start'].append(curr_utr_start)
                    data['utr_end'].append(curr_utr_end)
                    if int(curr_metadata_dict['exon_number']) == 1:
                        data['utr_type'].append(GenomicRegionTypes.FIVE_PRIME_UTR)
                    else:
                        data['utr_type'].append(GenomicRegionTypes.THREE_PRIME_UTR)
        self.df_utrs = pd.DataFrame(data)
        logger.info('Loaded %i UTRs in total.' % len(self.df_utrs))

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
                annotator=Annotators.GENCODE,
                annotator_version=self.version,
                region=GenomicRegionTypes.INTERGENIC,
                species=self.species
            )
            variant_call_annotations.append(variant_call_annotation)
        else:
            for _, row in df_genes_matched.iterrows():
                df_utrs_matched = self.df_utrs[
                    (self.df_utrs['gene_id'] == row['gene_id']) &
                    (self.df_utrs['utr_start'] <= position) &
                    (self.df_utrs['utr_end'] >= position)
                ]
                if len(df_utrs_matched) > 0:
                    region = df_utrs_matched['utr_type'].values[0]
                else:
                    df_exons_matched = self.df_exons[
                        (self.df_exons['gene_id'] == row['gene_id']) &
                        (self.df_exons['start'] <= position) &
                        (self.df_exons['end'] >= position)
                    ]
                    if len(df_exons_matched) > 0:
                        region = GenomicRegionTypes.EXONIC
                    else:
                        region = GenomicRegionTypes.INTRONIC
                variant_call_annotation = VariantCallAnnotation(
                    annotator=Annotators.GENCODE,
                    annotator_version=self.version,
                    gene_id=row['gene_id'],
                    gene_id_stable=row['gene_id_stable'],
                    gene_name=row['name'],
                    gene_strand=row['strand'],
                    gene_type=row['type'],
                    gene_version=row['version'],
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

