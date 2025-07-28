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


import gzip
import pandas as pd
import multiprocessing as mp
from collections import defaultdict
from dataclasses import dataclass, field, fields
from typing import List, Tuple
from .annotator import Annotator
from .constants import *
from .logging import get_logger
from .position_annotation import PositionAnnotation


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
        self._read_gtf_file_genes()        # add genes
        self._read_gtf_file_transcripts()  # add transcripts
        self._read_gtf_file_exons()        # add exons
        self._read_gtf_file_start_codons() # update start codon start and end positions
        self._read_gtf_file_stop_codons()  # update stop codon start and end positions
        self._read_gtf_file_utr()          # update UTR start and end positions

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

    def _read_gtf_file_genes(self):
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
        open_func = gzip.open if self.gtf_file.endswith('.gz') else open
        with open_func(self.gtf_file, 'rt') as f:
            for line in f:
                line = line.strip()
                if line.startswith('##'):
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
                    curr_level = int(curr_metadata_dict['level'])

                    if len(self.types) > 0 and curr_gene_type not in self.types:
                        continue
                    if len(self.levels) > 0 and curr_level not in self.levels:
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
                    data['level'].append(curr_level)
                    data['version'].append(curr_gene_version)
            self.df_genes = pd.DataFrame(data)
        self.__gene_ids = set()
        for gene_id in self.df_genes['gene_id'].values.tolist():
            self.__gene_ids.add(gene_id)
        logger.info('Loaded %i genes in total.' % len(self.df_genes))

    def _read_gtf_file_transcripts(self):
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
        open_func = gzip.open if self.gtf_file.endswith('.gz') else open
        with open_func(self.gtf_file, 'rt') as f:
            for line in f:
                line = line.strip()
                if line.startswith('##'):
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
                        curr_level = int(curr_metadata_dict['level'][0])
                    except:
                        curr_level = ''
                    try:
                        curr_transcript_support_level = int(curr_metadata_dict['transcript_support_level'][0])
                    except:
                        curr_transcript_support_level = ''
                    if curr_gene_id not in self.__gene_ids:
                        continue
                    if len(self.types) > 0 and curr_transcript_type not in self.types:
                        continue
                    if len(self.levels) > 0 and curr_level not in self.levels:
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
                    data['level'].append(curr_level)
                    data['support_level'].append(curr_transcript_support_level)
                    data['tags'].append(';'.join(curr_transcript_tags))
            self.df_transcripts = pd.DataFrame(data)
        self.__transcript_ids = set()
        for transcript_id in self.df_transcripts['transcript_id'].values.tolist():
            self.__transcript_ids.add(transcript_id)
        logger.info('Loaded %i transcripts in total.' % len(self.df_transcripts))

    def _read_gtf_file_exons(self):
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
            'level': [],
            'strand': [],
            'version': [],
            'tags': []
        }
        open_func = gzip.open if self.gtf_file.endswith('.gz') else open
        with open_func(self.gtf_file, 'rt') as f:
            for line in f:
                line = line.strip()
                if line.startswith('##'):
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
                    try:
                        curr_level = int(curr_metadata_dict['level'][0])
                    except:
                        curr_level = ''
                    if curr_gene_id not in self.__gene_ids:
                        continue
                    if curr_transcript_id not in self.__transcript_ids:
                        continue
                    if len(self.levels) > 0 and curr_level not in self.levels:
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
                    data['level'].append(curr_level)
                    data['strand'].append(curr_exon_strand)
                    data['version'].append(curr_exon_version)
                    data['tags'].append(';'.join(curr_exon_tags))
        self.df_exons = pd.DataFrame(data)
        logger.info('Loaded %i exons in total.' % len(self.df_exons))

    def _read_gtf_file_start_codons(self):
        """
        Read GENCODE GTF file and write rows to self.df_start_codons.
        """
        data = {
            'gene_id': [],
            'transcript_id': [],
            'start_codon_start': [],
            'start_codon_end': [],
            'level': []
        }
        open_func = gzip.open if self.gtf_file.endswith('.gz') else open
        with open_func(self.gtf_file, 'rt') as f:
            for line in f:
                line = line.strip()
                if line.startswith('##'):
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
                    try:
                        curr_level = int(curr_metadata_dict['level'][0])
                    except:
                        curr_level = ''
                    if len(self.levels) > 0 and curr_level not in self.levels:
                        continue
                    data['gene_id'].append(curr_gene_id)
                    data['transcript_id'].append(curr_transcript_id)
                    data['start_codon_start'].append(curr_start_codon_start)
                    data['start_codon_end'].append(curr_start_codon_end)
                    data['level'].append(curr_level)
        self.df_start_codons = pd.DataFrame(data)
        logger.info('Loaded %i start codons in total.' % len(self.df_start_codons))

    def _read_gtf_file_stop_codons(self):
        """
        Read GENCODE GTF file and write rows to self.df_stop_codons.
        """
        data = {
            'gene_id': [],
            'transcript_id': [],
            'stop_codon_start': [],
            'stop_codon_end': [],
            'level': []
        }
        open_func = gzip.open if self.gtf_file.endswith('.gz') else open
        with open_func(self.gtf_file, 'rt') as f:
            for line in f:
                line = line.strip()
                if line.startswith('##'):
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
                    try:
                        curr_level = int(curr_metadata_dict['level'][0])
                    except:
                        curr_level = ''
                    if len(self.levels) > 0 and curr_level not in self.levels:
                        continue
                    data['gene_id'].append(curr_gene_id)
                    data['transcript_id'].append(curr_transcript_id)
                    data['stop_codon_start'].append(curr_stop_codon_start)
                    data['stop_codon_end'].append(curr_stop_codon_end)
                    data['level'].append(curr_level)
        self.df_stop_codons = pd.DataFrame(data)
        logger.info('Loaded %i stop codons in total.' % len(self.df_stop_codons))

    def _read_gtf_file_utr(self):
        """
        Read GENCODE GTF file and write rows to self.df_utrs.
        """
        data = {
            'gene_id': [],
            'transcript_id': [],
            'utr_start': [],
            'utr_end': [],
            'utr_type': [],
            'level': []
        }
        open_func = gzip.open if self.gtf_file.endswith('.gz') else open
        with open_func(self.gtf_file, 'rt') as f:
            for line in f:
                line = line.strip()
                if line.startswith('##'):
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
                    try:
                        curr_level = int(curr_metadata_dict['level'][0])
                    except:
                        curr_level = ''
                    if len(self.levels) > 0 and curr_level not in self.levels:
                        continue
                    data['gene_id'].append(curr_gene_id)
                    data['transcript_id'].append(curr_transcript_id)
                    data['utr_start'].append(curr_utr_start)
                    data['utr_end'].append(curr_utr_end)
                    data['level'].append(curr_level)
                    if int(curr_metadata_dict['exon_number']) == 1:
                        data['utr_type'].append(GenomicRegionTypes.FIVE_PRIME_UTR)
                    else:
                        data['utr_type'].append(GenomicRegionTypes.THREE_PRIME_UTR)
        self.df_utrs = pd.DataFrame(data)
        logger.info('Loaded %i UTRs in total.' % len(self.df_utrs))

    def annotate_position(
            self,
            chromosome: str,
            position: int,
            strand: Strand
    ) -> List[PositionAnnotation]:
        """
        Annotate a position and return a list of PositionAnnotation objects.

        Parameters:
            chromosome              :   Chromosome.
            position                :   Position.

        Returns:
            List[PositionAnnotation]
        """
        if strand == Strand.FORWARD or strand == Strand.REVERSE:
            df_genes_matched = self.df_genes[
                (self.df_genes['chromosome'] == chromosome) &
                (self.df_genes['start'] <= position) &
                (self.df_genes['end'] >= position) &
                (self.df_genes['strand'] == strand)
            ]
        elif strand == Strand.BOTH or strand == Strand.UNKNOWN or pd.isna(strand):
            df_genes_matched = self.df_genes[
                (self.df_genes['chromosome'] == chromosome) &
                (self.df_genes['start'] <= position) &
                (self.df_genes['end'] >= position)
            ]
        else:
            raise Exception('Unknown strand: %s' % strand)

        if len(df_genes_matched) == 0:
            position_annotation = PositionAnnotation(
                annotator=Annotator.GENCODE,
                annotator_version=self.version,
                region=GenomicRegionTypes.INTERGENIC,
                species=self.species
            )
            return [position_annotation]
        position_annotations = []
        for _, row in df_genes_matched.iterrows():
            gene_id = row['gene_id']
            gene_id_stable = row['gene_id_stable']
            gene_name = row['name']
            gene_strand = row['strand']
            gene_type = row['type']
            gene_version = row['version']
            df_transcripts_matched = self.df_transcripts[
                (self.df_transcripts['gene_id'] == row['gene_id'])
            ]
            for _, row2 in df_transcripts_matched.iterrows():
                transcript_id = row2['transcript_id']
                transcript_id_stable = row2['transcript_id_stable']
                transcript_name = row2['name']
                transcript_strand = row2['strand']
                transcript_type = row2['type']
                transcript_version = row2['version']
                transcript_start = row2['start']
                transcript_end = row2['end']
                if transcript_start > position or transcript_end < position:
                    continue
                df_utrs_matched = self.df_utrs[
                    (self.df_utrs['gene_id'] == gene_id) &
                    (self.df_utrs['transcript_id'] == transcript_id) &
                    (self.df_utrs['utr_start'] <= position) &
                    (self.df_utrs['utr_end'] >= position)
                ]
                if len(df_utrs_matched) == 0:
                    df_exons_matched = self.df_exons[
                        (self.df_exons['gene_id'] == row['gene_id']) &
                        (self.df_exons['transcript_id'] == transcript_id) &
                        (self.df_exons['start'] <= position) &
                        (self.df_exons['end'] >= position)
                    ]
                    exon_id = ''
                    exon_id_stable = ''
                    if len(df_exons_matched) == 0:
                        region = GenomicRegionTypes.INTRONIC
                    else:
                        assert len(df_exons_matched) == 1
                        exon_id = df_exons_matched['exon_id'].values[0]
                        exon_id_stable = df_exons_matched['exon_id_stable'].values[0]
                        region = GenomicRegionTypes.EXONIC
                    position_annotation = PositionAnnotation(
                        annotator=Annotator.GENCODE,
                        annotator_version=self.version,
                        gene_id=gene_id,
                        gene_id_stable=gene_id_stable,
                        gene_name=gene_name,
                        gene_strand=gene_strand,
                        gene_type=gene_type,
                        gene_version=gene_version,
                        transcript_id=transcript_id,
                        transcript_id_stable=transcript_id_stable,
                        transcript_name=transcript_name,
                        transcript_strand=transcript_strand,
                        transcript_type=transcript_type,
                        transcript_version=transcript_version,
                        exon_id=exon_id,
                        exon_id_stable=exon_id_stable,
                        region=region,
                        species=self.species
                    )
                    position_annotations.append(position_annotation)
                else:
                    assert len(df_utrs_matched) == 1
                    region = df_utrs_matched['utr_type'].values[0]
                    position_annotation = PositionAnnotation(
                        annotator=Annotator.GENCODE,
                        annotator_version=self.version,
                        gene_id=gene_id,
                        gene_id_stable=gene_id_stable,
                        gene_name=gene_name,
                        gene_strand=gene_strand,
                        gene_type=gene_type,
                        gene_version=gene_version,
                        transcript_id=transcript_id,
                        transcript_id_stable=transcript_id_stable,
                        transcript_name=transcript_name,
                        transcript_strand=transcript_strand,
                        transcript_type=transcript_type,
                        transcript_version=transcript_version,
                        region=region,
                        species=self.species
                    )
                    position_annotations.append(position_annotation)
        if len(position_annotations) == 0:
            position_annotation = PositionAnnotation(
                annotator=Annotator.GENCODE,
                annotator_version=self.version,
                region=GenomicRegionTypes.INTERGENIC,
                species=self.species
            )
            position_annotations.append(position_annotation)
        return position_annotations

    def annotate_row(self, row) -> Tuple[str, List[PositionAnnotation], List[PositionAnnotation]]:
        # Position 1
        position_annotations_1 = self.annotate_position(
            chromosome=row['chromosome_1'],
            position=row['position_1'],
            strand=row['strand_1']
        )

        # Position 2
        position_annotations_2 = self.annotate_position(
            chromosome=row['chromosome_2'],
            position=row['position_2'],
            strand=row['strand_2']
        )

        return str(row['id']), position_annotations_1, position_annotations_2

    def annotate(self, df_variants: pd.DataFrame, num_processes: int) -> pd.DataFrame:
        """
        Annotate a Pandas DataFrame of standardized variants.

        Args:
            df_variants             :   Pandas DataFrame with the following columns:
                                        'id',
                                        'chromosome_1',
                                        'position_1',
                                        'strand_1',
                                        'chromosome_2',
                                        'position_2',
                                        'strand_2'
            num_processes           :   Number of processes.

        Returns:
            df_variants_annotated   :   Pandas DataFrame with the following columns appended:
                                        'annotator',
                                        'annotator_version',
                                        'gene_id',
                                        'gene_id_stable',
                                        'gene_name',
                                        'gene_strand',
                                        'gene_type',
                                        'gene_version',
                                        'transcript_id',
                                        'transcript_id_stable',
                                        'transcript_name',
                                        'transcript_strand',
                                        'transcript_type',
                                        'transcript_version',
                                        'exon_id',
                                        'exon_id_stable',
                                        'region',
                                        'species'
        """
        assert 'chromosome_1' in df_variants.columns
        assert 'chromosome_2' in df_variants.columns
        assert 'position_1' in df_variants.columns
        assert 'position_2' in df_variants.columns

        df_variants['id'] = df_variants['id'].astype(str)

        # Parallel annotation
        with mp.Pool(num_processes) as pool:
            results = pool.map(self.annotate_row, [row for _, row in df_variants.iterrows()])

        # Initialize annotation fields
        data = { 'id': [] }
        for f in fields(PositionAnnotation):
            data['position_1_%s' % f.name] = []
            data['position_2_%s' % f.name] = []

        # Populate results
        for variant_id, position_annotations_1, position_annotations_2 in results:
            position_annotations_1_dict = { key: [] for key in position_annotations_1[0].to_dict().keys() }
            position_annotations_2_dict = { key: [] for key in position_annotations_1[0].to_dict().keys() }
            for position_annotation in position_annotations_1:
                for key, value in position_annotation.to_dict().items():
                    position_annotations_1_dict[key].append(value)
            for position_annotation in position_annotations_2:
                for key, value in position_annotation.to_dict().items():
                    position_annotations_2_dict[key].append(value)
            position_annotations_1_data = { 'position_1_%s' % key: ';'.join(value) for key, value in position_annotations_1_dict.items()}
            position_annotations_2_data = { 'position_2_%s' % key: ';'.join(value) for key, value in position_annotations_2_dict.items()}
            data['id'].append(variant_id)
            for key, value in {**position_annotations_1_data, **position_annotations_2_data}.items():
                data[key].append(value)
        df_annotations = pd.DataFrame(data)

        df_variants_annotated = df_variants.merge(df_annotations, on='id', how='left')

        return df_variants_annotated