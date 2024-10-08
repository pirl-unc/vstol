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
The purpose of this python3 script is to implement the Ensembl dataclass.
"""


import pyensembl
from dataclasses import dataclass
from typing import List, Tuple
from .annotator import Annotator
from .constants import *
from .logging import get_logger
from .variant_call_annotation import VariantCallAnnotation
from .variant_call import VariantCall
from .variants_list import VariantsList


logger = get_logger(__name__)


@dataclass
class Ensembl(Annotator):
    release: int
    species: str
    _ensembl = None

    @property
    def ensembl(self) -> pyensembl.EnsemblRelease:
        if self._ensembl is None:
            self._ensembl = pyensembl.EnsemblRelease(release=self.release, species=self.species)
        return self._ensembl

    @property
    def source(self) -> str:
        return Annotators.ENSEMBL

    def annotate(self, variants_list: VariantsList, num_processes: int = 1) -> VariantsList:
        """
        Annotate a VariantsList object.

        Parameters:
            variants_list   :   VariantsList.

        Returns:
            VariantsList
        """
        for i in range(0, len(variants_list.variants)):
            for j in range(0, len(variants_list.variants[i].variant_calls)):
                position_1_annotations, position_2_annotations = self.annotate_variant_call_using_pyensembl(
                    variants_list.variants[i].variant_calls[j]
                )
                for annotation in position_1_annotations:
                    variants_list.variants[i].variant_calls[j].position_1_annotations.append(annotation)
                for annotation in position_2_annotations:
                    variants_list.variants[i].variant_calls[j].position_2_annotations.append(annotation)
        return variants_list

    def annotate_position_using_pyensembl(
            self,
            chromosome: str,
            position: int
    ) -> List[VariantCallAnnotation]:
        """
        Annotate a position using pyensembl and return a list of
        VariantAnnotation objects.

        Parameters:
            chromosome              :   Chromosome.
            position                :   Position.

        Returns:
            List[VariantCallAnnotation]
        """
        variant_call_annotations = []
        chromosome = chromosome.replace('chr', '')
        genes = self.ensembl.genes_at_locus(contig=chromosome, position=position)
        if len(genes) == 0:
            variant_call_annotation = VariantCallAnnotation(
                annotator=Annotators.ENSEMBL,
                region=GenomicRegionTypes.INTERGENIC,
                species=self.species,
                annotator_version=str(self.release)
            )
            variant_call_annotations.append(variant_call_annotation)
        else:
            for gene in genes:
                for transcript_id in self.ensembl.transcript_ids_of_gene_id(gene.gene_id):
                    transcript = self.ensembl.transcript_by_id(transcript_id)
                    if transcript.start > position or transcript.end < position:
                        continue
                    if transcript.contains_start_codon and transcript.contains_stop_codon: # protein-coding transcript
                        if transcript.strand == '+':
                            cds_start = transcript.start_codon_positions[0]
                            cds_end = transcript.stop_codon_positions[-1]
                            if position < cds_start or position > cds_end:
                                region = GenomicRegionTypes.UNTRANSLATED_REGION
                            else:
                                region = GenomicRegionTypes.INTRONIC
                        else:
                            cds_start = transcript.start_codon_positions[-1]
                            cds_end = transcript.stop_codon_positions[0]
                            if position > cds_start or position < cds_end:
                                region = GenomicRegionTypes.UNTRANSLATED_REGION
                            else:
                                region = GenomicRegionTypes.INTRONIC
                    else:
                        region = GenomicRegionTypes.INTRONIC
                    if region == GenomicRegionTypes.UNTRANSLATED_REGION:
                        variant_call_annotation = VariantCallAnnotation(
                            annotator=Annotators.ENSEMBL,
                            region=region,
                            species=self.species,
                            annotator_version=str(self.release),
                            gene_id=gene.gene_id,
                            gene_id_stable=gene.gene_id,
                            gene_name=gene.gene_name,
                            gene_strand=gene.strand,
                            gene_type=gene.biotype,
                            gene_version='',
                            transcript_id=transcript.transcript_id,
                            transcript_id_stable=transcript.transcript_id,
                            transcript_name=transcript.transcript_name,
                            transcript_strand=transcript.strand,
                            transcript_type=transcript.biotype
                        )
                    else:
                        exon_id = ''
                        for exon in transcript.exons:
                            if exon.start <= position <= exon.end:
                                region = GenomicRegionTypes.EXONIC
                                exon_id = exon.exon_id
                                break
                        variant_call_annotation = VariantCallAnnotation(
                            annotator=Annotators.ENSEMBL,
                            region=region,
                            species=self.species,
                            annotator_version=str(self.release),
                            gene_id=gene.gene_id,
                            gene_id_stable=gene.gene_id,
                            gene_name=gene.gene_name,
                            gene_strand=gene.strand,
                            gene_type=gene.biotype,
                            gene_version='',
                            transcript_id=transcript.transcript_id,
                            transcript_id_stable=transcript.transcript_id,
                            transcript_name=transcript.transcript_name,
                            transcript_strand=transcript.strand,
                            transcript_type=transcript.biotype,
                            transcript_version='',
                            exon_id=exon_id,
                            exon_id_stable=exon_id
                        )
                    variant_call_annotations.append(variant_call_annotation)
        return variant_call_annotations

    def annotate_variant_call_using_pyensembl(
            self,
            variant_call: VariantCall
    ) -> Tuple[List[VariantCallAnnotation], List[VariantCallAnnotation]]:
        """
        Annotate a VariantCall object and
        return two lists of VariantCallAnnotation objects.

        Parameters:
            variant_call            :   VariantCall.

        Returns:
            Tuple[position_1_annotations,position_2_annotations]
        """
        position_1_annotations = self.annotate_position_using_pyensembl(
            chromosome=variant_call.chromosome_1,
            position=variant_call.position_1
        )
        position_2_annotations = self.annotate_position_using_pyensembl(
            chromosome=variant_call.chromosome_2,
            position=variant_call.position_2
        )
        return position_1_annotations, position_2_annotations
