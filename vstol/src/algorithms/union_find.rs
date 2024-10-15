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


use std::collections::HashMap;


pub struct UnionFind {
    // key      =   child ID
    // value    =   parent ID
    pub parents: HashMap<String, String>,

    // key      =   parent ID
    // value    =   number of children
    sizes: HashMap<String, u32>
}

impl UnionFind {
    pub fn new() -> Self {
        UnionFind {
            parents: HashMap::new(),
            sizes: HashMap::new()
        }
    }

    pub fn get_clusters(&mut self) -> Vec<Vec<String>> {
        let keys: Vec<String> = self.parents.keys().cloned().collect();
        let mut map: HashMap<String, Vec<String>> = HashMap::new();
        for key in keys {
            let parent = self.find(&key);
            map.entry(parent)
                .or_insert_with(Vec::new)
                .push(key);
        }
        map.into_values().collect()
    }

    /// Get the size of a root ID.
    ///
    /// # Returns
    ///
    /// * Size of a root ID.
    fn get_size(&self, x: &str) -> u32 {
        *self.sizes.get(x).unwrap_or(&0)
    }

    /// Find the root ID of an element.
    ///
    /// # Returns
    ///
    /// * Root ID of the input element.
    pub fn find(&mut self, x: &str) -> String {
        if self.parents.contains_key(x) == false {
            self.parents.insert(x.to_string(), x.to_string());
            self.sizes.insert(x.to_string(), 1);
            return x.to_string();
        }

        // Compress path
        let mut path = Vec::new();
        let mut current = x;
        while let Some(parent) = self.parents.get(current) {
            if parent == current {
                break;
            }
            path.push(current.to_string());
            current = parent;
        }

        let root = current.to_string();
        for ancestor in path {
            self.parents.insert(ancestor, root.clone());
        }
        return root;
    }

    pub fn union(&mut self, x: &str, y: &str) {
        let root_x = self.find(x);
        let root_y = self.find(y);
        if root_x != root_y {
            let size_x = self.get_size(&root_x);
            let size_y = self.get_size(&root_y);
            if size_x < size_y {
                self.parents.insert(root_x.to_string(), root_y.to_string());
                self.sizes.entry(root_y.to_string()).and_modify(|e| *e += size_x);
            } else {
                self.parents.insert(root_y.clone(), root_x.clone());
                self.sizes.entry(root_x.clone()).and_modify(|e| *e += size_y);
            }
        }
    }
}