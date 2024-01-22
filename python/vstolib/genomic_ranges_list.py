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
The purpose of this python3 script is to implement the GenomicRangesList dataclass.
"""


import pandas as pd
from collections import defaultdict
from dataclasses import dataclass, field
from typing import List, Tuple, Dict
from .genomic_range import GenomicRange
from .logging import get_logger


logger = get_logger(__name__)


@dataclass(frozen=True)
class GenomicRangesList:
    # key   =   chromosome
    # value =   list of GenomicRange objects
    genomic_ranges_map = defaultdict(list)

    # key   =   GenomicRange.id
    # value =   (chromosome, index)
    _genomic_ranges_map: Dict[str, Tuple[str, int]] = field(default_factory=dict)

    def __post_init__(self):
        for key, val in self.genomic_ranges_map.items():
            i = 0
            for genomic_range in val:
                self._genomic_ranges_map[genomic_range.id] = (key, i)
                i += 1

    @property
    def size(self):
        size = 0
        for _ in self.genomic_ranges_map.values():
            size += 1
        return size

    @property
    def num_genomic_regions(self):
        return len(self.genomic_ranges_map.values())

    def add_genomic_range(self, genomic_range: GenomicRange):
        """
        Add a GenomicRange object.
        """
        curr_len = len(self.genomic_ranges_map[genomic_range.chromosome])
        self.genomic_ranges_map[genomic_range.chromosome].append(genomic_range)
        self._genomic_ranges_map[genomic_range.id] = (genomic_range.chromosome, curr_len)

    def find_overlaps(self, chromosome: str, start: int, end: int) -> List[GenomicRange]:
        """
        Find GenomicRange objects that overlap with query position.

        Parameters:
            chromosome      :   Chromosome.
            start           :   Start position.
            end             :   End position.

        Returns:
            List[GenomicRange]
        """
        # Get GenomicRange objects that match the query position
        genomic_ranges = []
        for genomic_range in self.genomic_ranges_map[chromosome]:
            if genomic_range.overlaps(chromosome=chromosome, start=start, end=end):
                genomic_ranges.append(genomic_range)
        return genomic_ranges

    def get_genomic_range(self, id: str) -> GenomicRange:
        chromosome = self._genomic_ranges_map[id][0]
        idx = self._genomic_ranges_map[id][1]
        return self.genomic_ranges_map[chromosome][idx]

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self.to_dataframe_row())

    def to_dataframe_row(self) -> Dict:
        data = defaultdict(list)
        for _, genomic_ranges  in self.genomic_ranges_map.items():
            for genomic_range in genomic_ranges:
                for key, values in genomic_range.to_dataframe_row().items():
                    for value in values:
                        data[key].append(value)
        return data

    def to_dict(self) -> Dict:
        data = {
            'genomic_ranges_map': {}
        }
        for chromosome, genomic_ranges in self.genomic_ranges_map.items():
            data['genomic_ranges_map'][chromosome] = []
            for genomic_range in genomic_ranges:
                data['genomic_ranges_map'][chromosome].append(genomic_range.to_dict())
        return data

    @staticmethod
    def load_dataframe(df: pd.DataFrame) -> 'GenomicRangesList':
        """
        Load a Pandas DataFrame and return a GenomicRangesList object.

        Parameters:
            df      :   Pandas DataFrame.

        Returns:
            GenomicRangesList
        """
        genomic_ranges_list = GenomicRangesList()
        for index, row in df.iterrows():
            chromosome = str(row['chromosome'])
            start = int(row['start'])
            end = int(row['end'])
            genomic_range = GenomicRange(
                chromosome=chromosome,
                start=start,
                end=end
            )
            genomic_ranges_list.add_genomic_range(genomic_range=genomic_range)
        return genomic_ranges_list

    @staticmethod
    def read_tsv_file(tsv_file: str,
                      low_memory: bool = True,
                      memory_map: bool = False) -> 'GenomicRangesList':
        """
        Read a TSV file and return a GenomicRangesList object.

        Parameters:
            tsv_file    :   TSV file.
            low_memory  :   Low memory (default: True).
            memory_map  :   Map memory (default: False).

        Returns:
            GenomicRangesList
        """
        df = pd.read_csv(tsv_file,
                         sep='\t',
                         low_memory=low_memory,
                         memory_map=memory_map)
        return GenomicRangesList.load_dataframe(df=df)
