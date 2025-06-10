extern crate nalgebra as na;
extern crate rand;

use na::{UnitVector3, Vector3};
use rand::{Rng, rng};

pub fn safe_normalize(vec: Vector3<f32>) -> UnitVector3<f32> {
    if vec.norm_squared() < 1e-10 {
        return UnitVector3::new_normalize(Vector3::x());
    }
    return UnitVector3::new_normalize(vec);
}

pub fn rand_vector3(min_val: f32, max_val: f32) -> Vector3<f32> {
    let mut rng = rng();

    return Vector3::new(
        rng.random_range(min_val..max_val),
        rng.random_range(min_val..max_val),
        rng.random_range(min_val..max_val),
    );
}

pub struct UnionFind {
    pub parent: Vec<usize>,
    rank: Vec<usize>,
}

impl UnionFind {
    pub fn new(size: usize) -> Self {
        Self { parent: (0..size).collect(), rank: vec![0; size] }
    }

    pub fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]);
        }
        self.parent[x]
    }

    pub fn union(&mut self, x: usize, y: usize) {
        let mut root_x = self.find(x);
        let mut root_y = self.find(y);

        if root_x == root_y {
            return;
        }

        if self.rank[root_x] < self.rank[root_y] {
            std::mem::swap(&mut root_x, &mut root_y);
        }

        self.parent[root_y] = root_x;
        if self.rank[root_x] == self.rank[root_y] {
            self.rank[root_x] += 1;
        }
    }
}
