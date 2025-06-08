extern crate nalgebra as na;
use crate::{
    rigid_body::AABB,
    world::{ColliderMap, OBBMap},
};
use core::f32;
use na::{UnitVector3, Vector3};
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
