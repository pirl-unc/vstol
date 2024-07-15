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
The purpose of this python3 script is to implement functions related to calculating metrics.
"""


import pysam


def calculate_average_alignment_score(
        bam_file: str,
        chromosome: str,
        position: int,
        window: int
) -> float:
    """
    Calculates the average alignment score in a genomic position.

    Args:
        bam_file        :   BAM file.
        chromosome      :   Chromosome.
        position        :   Position.
        window          :   Window (will be applied both upstream and downstream).

    Returns:
        Average alignment score.
    """
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    chromosome_length = bamfile.get_reference_length(chromosome)
    start = position - window
    end = position + window
    if start < 0:
        start = 0
    if end > chromosome_length:
        end = chromosome_length
    total = 0
    count = 0
    for read in bamfile.fetch(chromosome, start, end):
        total += read.mapping_quality
        count += 1
    if count == 0:
        return -1.0
    return float(total) / float(count)

