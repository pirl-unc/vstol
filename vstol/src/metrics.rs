// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::collections::{HashMap};


pub fn calculate_average_alignment_scores(
    bam_file: &str,
    regions: &Vec<(String, u32, u32)>,
    num_threads: usize,
) -> HashMap<(String, u32, u32), f64> {
    // Step 1. Read the BAM file
    let reader = bam::IndexedReader::from_path(bam_file).unwrap();
    let header = reader.header();

    // Step 2. Get the regions with chromosome IDs
    let mut regions_reformatted: Vec<(u32, u32, u32)> = Vec::new();
    let mut chromosomes_map: HashMap<u32, String> = HashMap::new();
    for (chromosome, start, end) in regions.iter() {
        if let Some(chromosome_id) = header.reference_id(chromosome) {
            regions_reformatted.push((chromosome_id, *start, *end));
            chromosomes_map.insert(chromosome_id, chromosome.to_string());
        } else {
            panic!("The chromosome ID cannot be found for {}", chromosome);
        }
    }

    // Step 3. Split the regions into roughly equal-sized chunks
    let chunk_size = (regions_reformatted.len() + num_threads - 1) / num_threads;
    let regions_reformatted_chunks: Vec<_> = regions_reformatted.chunks(chunk_size).collect();

    // Step 4. Calculate the average alignment score for each position
    let thread_pool = ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let scores: Vec<(u32, u32, u32, f64)> = thread_pool.install(|| {
        regions_reformatted_chunks
            .par_iter()
            .flat_map(|chunk| {
                let mut local_reader = bam::IndexedReader::from_path(bam_file).unwrap();
                chunk
                    .iter()
                    .map(|(chromosome, start, end)| {
                        let mut total: u32 = 0;
                        let mut count: u32 = 0;
                        for record in local_reader.fetch(&bam::Region::new(*chromosome, *start, *end)).unwrap() {
                            total += record.unwrap().mapq() as u32;
                            count += 1;
                        }
                        let score = if count > 0 {
                            total as f64 / count as f64
                        } else {
                            -1.0
                        };
                        (*chromosome, *start, *end, score)
                    })
                    .collect::<Vec<(u32, u32, u32, f64)>>()
            })
            .collect()
    });

    // Step 5. Store the average alignment scores as a HashMap
    let mut scores_map: HashMap<(String, u32, u32), f64> = HashMap::new();
    for (chromosome_id, start, end, score) in scores.iter() {
        if let Some(chromosome_name) = chromosomes_map.get(chromosome_id) {
            scores_map.insert((chromosome_name.clone(), *start, *end), *score);
        }
    }

    scores_map
}
