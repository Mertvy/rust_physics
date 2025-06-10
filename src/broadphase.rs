extern crate nalgebra as na;
use crate::{
    rigid_body::AABB,
    world::{ColliderMap, OBBMap},
};
use core::f32;
use na::{UnitVector3, Vector3};
use std::collections::{HashMap, HashSet};
use std::usize;

pub type ArrayList<T> = Vec<T>;
pub struct ColliderPair(pub usize, pub usize);

pub trait BroadPhase {
    fn add(&mut self, aabb: AABB);
    fn update(&mut self, colliders: &ColliderMap, obbs: &OBBMap);
    fn query_potential_collisions(&self) -> ArrayList<ColliderPair>;
    fn query_point(&self, point: &Vector3<f32>) -> Option<usize>;
    fn query_aabb(&self, aabb: &AABB) -> ArrayList<usize>;
    fn query_ray(&self, ray: &Ray3, colliders: &ColliderMap) -> RayCastResult;
}

pub struct Ray3 {
    pub pos: Vector3<f32>,
    pub dir: UnitVector3<f32>,
}

pub struct RayCastResult {
    hit: bool,
    collider_id: usize,
    position: Vector3<f32>,
    normal: Vector3<f32>,
}

pub struct NSquared {
    pub aabb_list: ArrayList<AABB>,
}

impl BroadPhase for NSquared {
    fn add(&mut self, aabb: AABB) {
        self.aabb_list.push(aabb);
    }

    fn update(&mut self, colliders: &ColliderMap, obbs: &OBBMap) {
        for aabb in &mut self.aabb_list {
            let obb = obbs
                .get(&colliders.get(&aabb.collider_id).expect("AABB has not collider").obb_id)
                .expect("Collider has no OBB");
            aabb.update_from_obb(obb);
        }
    }

    fn query_potential_collisions(&self) -> ArrayList<ColliderPair> {
        let mut collider_pair_list = ArrayList::with_capacity(self.aabb_list.capacity() / 4);
        for i in 0..self.aabb_list.len() {
            for j in i + 1..self.aabb_list.len() {
                let aabb1 = self.aabb_list.get(i).expect("AABB list out of whack");
                let aabb2 = self.aabb_list.get(j).expect("AABB list out of whack");

                if aabb1.collides(aabb2) {
                    collider_pair_list.push(ColliderPair(aabb1.collider_id, aabb2.collider_id));
                }
            }
        }
        return collider_pair_list;
    }

    fn query_point(&self, point: &Vector3<f32>) -> Option<usize> {
        for aabb in &self.aabb_list[..] {
            if aabb.contains(point) {
                return Some(aabb.collider_id);
            }
        }
        return None;
    }

    fn query_aabb(&self, aabb: &AABB) -> ArrayList<usize> {
        let mut query_list = ArrayList::new();
        for other in self.aabb_list.iter() {
            if aabb.collides(other) {
                query_list.push(other.collider_id);
            }
        }
        return query_list;
    }

    fn query_ray(&self, ray: &Ray3, colliders: &ColliderMap) -> RayCastResult {
        let mut result = RayCastResult {
            hit: false,
            collider_id: usize::MAX,
            position: Vector3::zeros(),
            normal: Vector3::zeros(),
        };

        let mut candidates: ArrayList<usize> = ArrayList::with_capacity(self.aabb_list.len());

        for aabb in &self.aabb_list {
            if aabb.intersects_ray(ray) {
                candidates.push(aabb.collider_id);
            }
        }

        let mut min_t = f32::INFINITY;
        for collider_id in candidates {
            let collider = colliders.get(&collider_id).expect("Requested collider does not exist");
            let (hit, t, normal) = collider.intersects_ray(ray);
            if !hit || t > min_t {
                continue;
            };
            min_t = t;
            result.hit = true;
            result.collider_id = collider_id;
            result.position = ray.pos + ray.dir.scale(t);
            result.normal = normal;
        }
        return result;
    }
}

fn world_to_cell(pos: &Vector3<f32>, cell_size: f32) -> (i32, i32, i32) {
    (
        (pos.x / cell_size).floor() as i32,
        (pos.y / cell_size).floor() as i32,
        (pos.z / cell_size).floor() as i32,
    )
}

pub struct SpatialHashing {
    cell_size: f32,
    spatial_map: HashMap<(i32, i32, i32), Vec<usize>>,
}

impl SpatialHashing {
    pub fn new(cell_size: f32) -> Self {
        Self { cell_size, spatial_map: HashMap::new() }
    }
}

impl BroadPhase for SpatialHashing {
    fn update(&mut self, colliders: &ColliderMap, obbs: &OBBMap) {
        self.spatial_map.clear();

        for (collider_id, collider) in colliders {
            let aabb = AABB::new_from_obb(*collider_id, &obbs[&collider.obb_id]);
            let min = world_to_cell(&aabb.global_min, self.cell_size);
            let max = world_to_cell(&aabb.global_max, self.cell_size);

            for x in min.0..=max.0 {
                for y in min.1..=max.1 {
                    for z in min.2..=max.2 {
                        self.spatial_map.entry((x, y, z)).or_default().push(*collider_id);
                    }
                }
            }
        }
    }

    fn query_potential_collisions(&self) -> ArrayList<ColliderPair> {
        let mut pair_set = HashSet::new();
        let mut pairs = ArrayList::new();
        for cell_contents in self.spatial_map.values() {
            for i in 0..cell_contents.len() {
                for j in (i + 1)..cell_contents.len() {
                    let id1 = cell_contents[i];
                    let id2 = cell_contents[j];
                    let key = (id1.min(id2), id1.max(id2));
                    if pair_set.insert(key) {
                        pairs.push(ColliderPair(id1, id2));
                    }
                }
            }
        }

        return pairs;
    }

    fn add(&mut self, aabb: AABB) {}

    fn query_aabb(&self, aabb: &AABB) -> ArrayList<usize> {
        todo!();
    }

    fn query_point(&self, point: &Vector3<f32>) -> Option<usize> {
        todo!()
    }

    fn query_ray(&self, ray: &Ray3, colliders: &ColliderMap) -> RayCastResult {
        todo!()
    }
}
