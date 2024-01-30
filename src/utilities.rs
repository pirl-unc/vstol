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


use std::collections::{HashMap, HashSet};


pub fn calculate_max(values: Vec<f64>) -> f64 {
    let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    return max;
}

pub fn calculate_mean(values: Vec<f64>) -> f64 {
    let sum: f64 = values.iter().sum();
    let average = sum as f64 / values.len() as f64;
    return average;
}

pub fn calculate_median(values: Vec<f64>) -> f64 {
    let mut values_: Vec<f64> = values.clone();
    values_.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut median: f64 = 0.0;
    let len = values_.len();
    if len % 2 == 0 {
        // If the number of elements is even,
        // take the average of the two middle values
        let middle1 = len / 2 - 1;
        let middle2 = len / 2;
        median = (values_[middle1] + values_[middle2]) / 2.0;
    } else {
        // If the number of elements is odd,
        // take the middle value
        let middle = len / 2;
        median = values_[middle];
    }
    return median;
}

pub fn calculate_min(values: Vec<f64>) -> f64 {
    let min = values.iter().cloned().fold(f64::INFINITY, f64::min);
    return min;
}

pub fn find_clusters(pairs: Vec<(String, String)>) -> Vec<HashSet<String>> {
    let mut adjacency_map: HashMap<String, HashSet<String>> = HashMap::new();
    for (element_1, element_2) in pairs.iter() {
        adjacency_map
            .entry(element_1.clone())
            .or_insert_with(HashSet::new)
            .insert(element_2.clone());
        adjacency_map
            .entry(element_2.clone())
            .or_insert_with(HashSet::new)
            .insert(element_1.clone());
    }
    let mut visited: HashSet<String> = HashSet::new();
    let mut clusters: Vec<HashSet<String>> = Vec::new();
    fn dfs(element: String,
           cluster: &mut HashSet<String>,
           visited: &mut HashSet<String>,
           adjacency_map: &HashMap<String,
           HashSet<String>>) {
        visited.insert(element.clone());
        cluster.insert(element.clone());
        if let Some(neighbors) = adjacency_map.get(&element) {
            for neighbor in neighbors {
                if !visited.contains(neighbor) {
                    dfs(neighbor.clone(), cluster, visited, adjacency_map);
                }
            }
        }
    }
    for element in adjacency_map.keys() {
        if visited.contains(element) == false {
            let mut cluster = HashSet::new();
            dfs(element.clone(), &mut cluster, &mut visited, &adjacency_map);
            clusters.push(cluster);
        }
    }
    return clusters;
}
