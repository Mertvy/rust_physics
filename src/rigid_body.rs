extern crate nalgebra as na;
use crate::{
    broadphase::Ray3,
    utils::safe_normalize,
    world::{ColliderMap, RigidBodyMap},
};
use core::f32;
use na::{Matrix3, UnitQuaternion, Vector3};
use std::array;

type ArrayList<T> = Vec<T>;

#[derive(Debug, Clone)]
pub struct RigidBody {
    pub id: usize,
    pub mass: f32,
    pub inv_mass: f32,
    pub global_inv_moment_inertia: Matrix3<f32>,
    local_inv_moment_inertia: Matrix3<f32>,

    pub global_center_mass: Vector3<f32>,
    local_center_mass: Vector3<f32>,

    pub global_orientation: UnitQuaternion<f32>,
    pub lin_velocity: Vector3<f32>,
    pub ang_velocity: Vector3<f32>,

    force_accumulator: Vector3<f32>,
    torque_accumulator: Vector3<f32>,

    collider_id_list: ArrayList<usize>,

    pub is_static: bool,
}

impl RigidBody {
    pub fn new_body(
        id: usize,
        position: Vector3<f32>,
        orientation: UnitQuaternion<f32>,
        is_static: bool,
    ) -> RigidBody {
        return RigidBody {
            id: id,
            mass: 0.,
            inv_mass: f32::INFINITY,
            global_inv_moment_inertia: Matrix3::zeros(),
            local_inv_moment_inertia: Matrix3::zeros(),

            global_center_mass: position,
            local_center_mass: Vector3::zeros(),

            global_orientation: orientation,
            lin_velocity: Vector3::zeros(),
            ang_velocity: Vector3::zeros(),

            force_accumulator: Vector3::zeros(),
            torque_accumulator: Vector3::zeros(),

            collider_id_list: ArrayList::new(),

            is_static: is_static,
        };
    }

    pub fn update(&mut self, dt: f32) {
        if self.is_static {
            return;
        }
        self.lin_velocity *= 0.998;
        self.lin_velocity += (self.force_accumulator * self.inv_mass) * dt;
        self.force_accumulator = Vector3::zeros();

        self.ang_velocity *= 0.998;
        self.ang_velocity += (self.global_inv_moment_inertia * self.torque_accumulator) * dt;
        self.torque_accumulator = Vector3::zeros();

        self.global_center_mass += self.lin_velocity * dt;

        if self.ang_velocity.norm() > 1e-6 {
            let rot_axis = safe_normalize(self.ang_velocity);
            let rot_angle = self.ang_velocity.norm() * dt;
            let rot_operator = UnitQuaternion::from_axis_angle(&rot_axis, rot_angle);
            self.global_orientation = rot_operator * self.global_orientation;
            self.global_orientation.renormalize();
        }

        let rot_mat = self.global_orientation.to_rotation_matrix();
        self.global_inv_moment_inertia =
            rot_mat * self.local_inv_moment_inertia * rot_mat.transpose();
    }

    pub fn add_collider(&mut self, collider_id: usize, colliders: &ColliderMap) {
        self.collider_id_list.push(collider_id);

        self.local_center_mass = Vector3::zeros();
        self.mass = 0.;

        for id in &self.collider_id_list {
            let collider = colliders.get(id).expect("Collider does not exist");
            self.mass += collider.mass;
            self.local_center_mass += collider.local_center_mass * collider.mass;
        }

        self.inv_mass = 1. / self.mass;
        self.local_center_mass *= self.inv_mass;

        // Parallel Axis Theorem: https://en.wikipedia.org/wiki/Parallel_axis_theorem#Tensor_generalization

        let mut local_inertia_tensor: Matrix3<f32> = Matrix3::zeros();
        for id in &self.collider_id_list {
            let collider = colliders.get(id).expect("Collider does not exist");
            let r: Vector3<f32> = self.local_center_mass - collider.local_center_mass;
            let dist_squared = r.dot(&r);
            let outer_prod: Matrix3<f32> = r * r.transpose();
            local_inertia_tensor += collider.local_inertia_tensor
                + collider.mass * (dist_squared * Matrix3::identity() - outer_prod);
        }
        self.local_inv_moment_inertia =
            local_inertia_tensor.try_inverse().expect("Inertia tensor not invertible");
    }

    pub fn local_to_global_pos(&self, pos: &Vector3<f32>) -> Vector3<f32> {
        return self.global_orientation.transform_vector(pos) + self.global_center_mass;
    }

    pub fn global_to_local_pos(&self, pos: &Vector3<f32>) -> Vector3<f32> {
        let v: Vector3<f32> = pos - &self.global_center_mass;
        return self.global_orientation.inverse_transform_vector(&v);
    }

    pub fn local_to_global_vec(&self, vec: &Vector3<f32>) -> Vector3<f32> {
        return self.global_orientation.transform_vector(vec);
    }

    pub fn global_to_local_vec(&self, vec: &Vector3<f32>) -> Vector3<f32> {
        return self.global_orientation.inverse_transform_vector(vec);
    }

    pub fn apply_force_local(&mut self, force: &Vector3<f32>, appl_pos: &Vector3<f32>) {
        let global_force: Vector3<f32> = self.local_to_global_vec(force);
        let global_appl_pos: Vector3<f32> = self.local_to_global_pos(appl_pos);
        self.apply_force_global(&global_force, &global_appl_pos);
    }

    pub fn apply_force_global(&mut self, force: &Vector3<f32>, appl_pos: &Vector3<f32>) {
        self.force_accumulator += force;
        self.torque_accumulator += (appl_pos - self.global_center_mass).cross(force)
    }
}

#[derive(Debug)]
pub enum ColliderShape {
    Sphere { radius: f32 },
    Box { x_len: f32, y_len: f32, z_len: f32 },
}

pub struct Collider {
    pub id: usize,
    pub body_id: usize,

    pub obb_id: usize,

    pub shape: ColliderShape,
    pub mass: f32,
    local_inertia_tensor: Matrix3<f32>,

    local_center_mass: Vector3<f32>,
    local_orientation: UnitQuaternion<f32>,
}

impl Collider {
    pub fn new_collider(
        id: usize,
        body_id: usize,
        obb_id: usize,
        shape: ColliderShape,
        mass: f32,
        local_center_mass: Vector3<f32>,
        local_orientation: UnitQuaternion<f32>,
    ) -> Collider {
        let local_inertia_tensor = match shape {
            ColliderShape::Sphere { radius } => {
                (2. / 5.) * mass * radius * radius * Matrix3::identity()
            }
            ColliderShape::Box { x_len, y_len, z_len } => {
                (1. / 12.)
                    * mass
                    * Matrix3::from_diagonal(&Vector3::new(
                        y_len * y_len + z_len * z_len,
                        x_len * x_len + z_len * z_len,
                        x_len * x_len + y_len * y_len,
                    ))
            }
        };

        Collider {
            id: id,
            body_id: body_id,
            obb_id: obb_id,
            shape: shape,
            mass: mass,
            local_inertia_tensor: local_inertia_tensor,
            local_center_mass: local_center_mass,
            local_orientation: local_orientation,
        }
    }

    pub fn support(&self, local_dir: &Vector3<f32>) -> Vector3<f32> {
        match &self.shape {
            ColliderShape::Sphere { radius } => {
                return self.local_center_mass + *radius * local_dir.normalize();
            }
            ColliderShape::Box { x_len, y_len, z_len } => {
                return self.local_center_mass
                    + Vector3::new(
                        if local_dir.x >= 0. { x_len / 2. } else { -x_len / 2. },
                        if local_dir.y >= 0. { y_len / 2. } else { -y_len / 2. },
                        if local_dir.z >= 0. { z_len / 2. } else { -z_len / 2. },
                    );
            }
        };
    }

    pub fn intersects_ray(&self, ray: &Ray3) -> (bool, f32, Vector3<f32>) {
        return (false, f32::INFINITY, Vector3::zeros());
    }
}

fn is_overlapping(int1: (f32, f32), int2: (f32, f32)) -> bool {
    return !(int1.1 < int2.0 || int1.0 > int2.1);
}

pub struct OBB {
    pub id: usize,
    pub collider_id: usize,
    pub local_min: Vector3<f32>,
    pub local_max: Vector3<f32>,
    pub global_pos: Vector3<f32>,
    pub global_orientation: UnitQuaternion<f32>,
}

impl OBB {
    pub fn new(
        id: usize,
        collider_id: usize,
        colliders: &ColliderMap,
        rigid_bodies: &RigidBodyMap,
    ) -> OBB {
        let mut obb = OBB {
            id: id,
            collider_id: collider_id,
            local_min: Vector3::zeros(),
            local_max: Vector3::zeros(),
            global_pos: Vector3::zeros(),
            global_orientation: UnitQuaternion::identity(),
        };
        obb.update(colliders, rigid_bodies);
        return obb;
    }

    pub fn update(&mut self, colliders: &ColliderMap, rigid_bodies: &RigidBodyMap) {
        let collider = colliders.get(&self.collider_id).expect("OBB has no collider");
        let rigid_body = rigid_bodies.get(&collider.body_id).expect("Collider has no body");

        self.global_pos = rigid_body.local_to_global_pos(&collider.local_center_mass);
        self.global_orientation = rigid_body.global_orientation * collider.local_orientation;

        match &collider.shape {
            ColliderShape::Box { x_len, y_len, z_len } => {
                self.local_min = Vector3::new(-x_len / 2., -y_len / 2., -z_len / 2.);
                self.local_max = Vector3::new(x_len / 2., y_len / 2., z_len / 2.);
            }
            ColliderShape::Sphere { radius } => {
                self.local_min = Vector3::new(-radius, -radius, -radius);
                self.local_max = Vector3::new(*radius, *radius, *radius);
            }
        }
    }

    pub fn contains(&self, point: &Vector3<f32>) -> bool {
        let p = self.global_orientation.inverse_transform_vector(&(point - self.global_pos));
        return (self.local_min.x <= p.x && p.x <= self.local_max.x)
            && (self.local_min.y <= p.y && p.y <= self.local_max.y)
            && (self.local_min.z <= p.z && p.z <= self.local_max.z);
    }

    pub fn collides(&self, obb: &OBB) -> bool {
        // SEPARATING AXIS THEOREM https://dyn4j.org/2010/01/sat/

        let axes1 = self.compute_normal_axes();
        let axes2 = obb.compute_normal_axes();

        let global_vertices1 = self.get_global_vertices();
        let global_vertices2 = obb.get_global_vertices();

        for axis in axes1 {
            let int1 = OBB::project_onto_axis(&global_vertices1, &axis);
            let int2 = OBB::project_onto_axis(&global_vertices2, &axis);
            if !is_overlapping(int1, int2) {
                return false;
            }
        }

        for axis in axes2 {
            let int1 = OBB::project_onto_axis(&global_vertices1, &axis);
            let int2 = OBB::project_onto_axis(&global_vertices2, &axis);
            if !is_overlapping(int1, int2) {
                return false;
            }
        }

        for axis1 in axes1 {
            for axis2 in axes2 {
                let cross = axis1.cross(&axis2);
                if cross.norm_squared() < 1e-6 {
                    continue;
                }
                let int1 = OBB::project_onto_axis(&global_vertices1, &cross);
                let int2 = OBB::project_onto_axis(&global_vertices2, &cross);
                if !is_overlapping(int1, int2) {
                    return false;
                }
            }
        }

        return true;
    }

    fn compute_normal_axes(&self) -> [Vector3<f32>; 3] {
        let global_x = self.global_orientation.transform_vector(&Vector3::new(1., 0., 0.));
        let global_y = self.global_orientation.transform_vector(&Vector3::new(0., 1., 0.));
        let global_z = self.global_orientation.transform_vector(&Vector3::new(0., 0., 1.));
        return [global_x, global_y, global_z];
    }

    fn project_onto_axis(global_vertices: &[Vector3<f32>; 8], axis: &Vector3<f32>) -> (f32, f32) {
        let mut min_val = f32::INFINITY;
        let mut max_val = f32::NEG_INFINITY;

        for vertex in global_vertices {
            let dot = vertex.dot(axis);

            min_val = f32::min(min_val, dot);
            max_val = f32::max(max_val, dot);
        }

        return (min_val, max_val);
    }

    fn get_local_vertices(&self) -> [Vector3<f32>; 8] {
        return array::from_fn(|i| {
            let x = if i & 1 == 0 { self.local_min.x } else { self.local_max.x };
            let y = if i & 2 == 0 { self.local_min.y } else { self.local_max.y };
            let z = if i & 4 == 0 { self.local_min.z } else { self.local_max.z };
            return Vector3::new(x, y, z);
        });
    }

    pub fn get_global_vertices(&self) -> [Vector3<f32>; 8] {
        let local_axes = self.get_local_vertices();
        return array::from_fn(|i| {
            return self.global_pos + self.global_orientation.transform_vector(&local_axes[i]);
        });
    }
}

pub struct AABB {
    pub collider_id: usize,
    pub global_min: Vector3<f32>,
    pub global_max: Vector3<f32>,
}

impl AABB {
    pub fn new_from_obb(collider_id: usize, obb: &OBB) -> AABB {
        let mut aabb = AABB {
            collider_id: collider_id,
            global_min: Vector3::zeros(),
            global_max: Vector3::zeros(),
        };
        aabb.update_from_obb(obb);

        return aabb;
    }

    pub fn update_from_obb(&mut self, obb: &OBB) {
        let obb_vertices = obb.get_global_vertices();

        self.global_min = Vector3::new(f32::INFINITY, f32::INFINITY, f32::INFINITY);
        self.global_max = Vector3::new(f32::NEG_INFINITY, f32::NEG_INFINITY, f32::NEG_INFINITY);

        for vertex in &obb_vertices {
            self.global_min = self.global_min.inf(&vertex);
            self.global_max = self.global_max.sup(&vertex);
        }
    }

    pub fn contains(&self, point: &Vector3<f32>) -> bool {
        return (self.global_min.x <= point.x && point.x <= self.global_max.x)
            && (self.global_min.y <= point.y && point.y <= self.global_max.y)
            && (self.global_min.z <= point.z && point.z <= self.global_max.z);
    }

    pub fn intersects_ray(&self, ray: &Ray3) -> bool {
        todo!()
    }

    pub fn collides(&self, aabb: &AABB) -> bool {
        return is_overlapping(
            (self.global_min.x, self.global_max.x),
            (aabb.global_min.x, aabb.global_max.x),
        ) && is_overlapping(
            (self.global_min.y, self.global_max.y),
            (aabb.global_min.y, aabb.global_max.y),
        ) && is_overlapping(
            (self.global_min.z, self.global_max.z),
            (aabb.global_min.z, aabb.global_max.z),
        );
    }
}
