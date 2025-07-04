extern crate nalgebra as na;
extern crate rand;

use na::{UnitVector3, Vector3};
use rand::{Rng, rng};
use rayon::prelude::*;

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
    pub representatives: Vec<usize>,
    rank: Vec<usize>,
}

impl UnionFind {
    pub fn new(size: usize) -> Self {
        Self { representatives: (0..size).collect(), rank: vec![0; size] }
    }

    pub fn find(&mut self, x: usize) -> usize {
        if self.representatives[x] != x {
            self.representatives[x] = self.find(self.representatives[x]);
        }
        self.representatives[x]
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

        self.representatives[root_y] = root_x;
        if self.rank[root_x] == self.rank[root_y] {
            self.rank[root_x] += 1;
        }
    }
}

pub struct CustomMap<T> {
    storage: Vec<Option<T>>,
    size: usize,
}

impl<T> CustomMap<T> {
    pub fn new() -> Self {
        Self { storage: Vec::new(), size: 0 }
    }

    pub fn insert(&mut self, id: usize, value: T) {
        if id >= self.storage.len() {
            self.storage.resize_with(id + 1, || None);
        }
        self.storage[id] = Some(value);
        self.size += 1;
    }

    pub fn get(&self, id: &usize) -> Option<&T> {
        self.storage.get(*id)?.as_ref()
    }

    pub fn values(&self) -> impl Iterator<Item = &T> {
        self.storage.iter().filter_map(|x| x.as_ref())
    }

    pub fn get_mut(&mut self, id: &usize) -> Option<&mut T> {
        self.storage.get_mut(*id)?.as_mut()
    }

    pub fn values_mut(&mut self) -> impl Iterator<Item = &mut T> {
        self.storage.iter_mut().filter_map(|x| x.as_mut())
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut T> {
        self.storage.iter_mut().filter_map(|v| v.as_mut())
    }

    pub fn par_iter_mut(&mut self) -> impl ParallelIterator<Item = &mut T>
    where
        T: Send,
    {
        self.storage.par_iter_mut().filter_map(|v| v.as_mut())
    }

    pub fn len(&self) -> usize {
        return self.size;
    }
}

impl<T> Default for CustomMap<T> {
    fn default() -> Self {
        CustomMap { storage: Vec::new(), size: 0 }
    }
}
