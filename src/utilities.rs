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

struct UnionFind {
    parent: HashMap<String, String>
}

impl UnionFind {
    fn new() -> Self {
        UnionFind {
            parent: HashMap::new(),
        }
    }

    fn find(&mut self, x: &str) -> String {
        let parent = self.parent.entry(x.to_string()).or_insert_with(|| x.to_string()).clone();
        if parent != x {
            let root = self.find(&parent);
            self.parent.insert(x.to_string(), root.clone());
            root
        } else {
            x.to_string()
        }
    }

    fn union(&mut self, x: &str, y: &str) {
        let root_x = self.find(x);
        let root_y = self.find(y);
        if root_x != root_y {
            self.parent.insert(root_x, root_y);
        }
    }
}

pub fn find_clusters(pairs: HashSet<(String, String)>) -> Vec<Vec<String>> {
    let mut uf = UnionFind::new();
    for (x, y) in pairs {
        uf.union(&x, &y);
    }
    let ids: Vec<String> = uf.parent.keys().cloned().collect();
    let mut clusters_map: HashMap<String, Vec<String>> = HashMap::new();
    for id in ids {
        let root = uf.find(&id);
        clusters_map.entry(root).or_insert_with(Vec::new).push(id);
    }
    let mut clusters: Vec<Vec<String>> = Vec::new();
    for (root, ids) in clusters_map {
        clusters.push(ids);
    }
    return clusters;
}
