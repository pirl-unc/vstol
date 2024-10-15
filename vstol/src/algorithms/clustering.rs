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


use std::collections::{HashMap,HashSet};
use crate::algorithms::union_find::UnionFind;


/// Find clusters.
///
/// # Arguments
///
/// * `pairs` is a HashSet of ID pairs.
///
/// # Returns
///
/// * A vector of clusters of IDs.
pub fn find_clusters(pairs: HashSet<(String, String)>) -> Vec<Vec<String>> {
    let mut uf = UnionFind::new();
    for (x, y) in pairs {
        uf.union(&x, &y);
    }
    let ids: Vec<String> = uf.parents.keys().cloned().collect();
    let mut clusters_map: HashMap<String, Vec<String>> = HashMap::new();
    for id in ids {
        let root = uf.find(&id);
        clusters_map.entry(root).or_insert_with(Vec::new).push(id);
    }
    let mut clusters: Vec<Vec<String>> = Vec::new();
    for (_, ids) in clusters_map {
        clusters.push(ids);
    }
    clusters
}
