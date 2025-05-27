use crate::rigid_body::{Collider, RigidBody};
use std::collections::HashMap;

pub type RigidBodyMap = HashMap<usize, RigidBody>;
pub type ColliderMap = HashMap<usize, Collider>;

pub struct World {
    pub rigid_bodies: RigidBodyMap,
    pub colliders: ColliderMap,
}
