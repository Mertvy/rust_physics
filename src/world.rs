use nalgebra::{ArrayStorage, Matrix3, U12, UnitQuaternion, Vector, Vector3};

use crate::broadphase::{ArrayList, BroadPhase, ColliderPair, NSquared};
use crate::collision::{Contact, expanding_polytope, gjk};
use crate::constraint::VelocityConstraint;
use crate::rigid_body::{AABB, Collider, ColliderShape, OBB, RigidBody};
use core::f32;
use std::collections::HashMap;

pub type RigidBodyMap = HashMap<usize, RigidBody>;
pub type ColliderMap = HashMap<usize, Collider>;
pub type OBBMap = HashMap<usize, OBB>;

const dt: f32 = 1. / 60.;
const GRAVITY: Vector3<f32> = Vector3::new(0., -30., 0.);

pub struct World {
    pub rigid_bodies: RigidBodyMap,
    next_body_id: usize,

    pub colliders: ColliderMap,
    next_collider_id: usize,

    pub obbs: OBBMap,
    next_obb_id: usize,

    broadphase: NSquared,
}

impl World {
    pub fn step(&mut self) {
        self.apply_gravity();

        let contacts = self.generate_contacts();

        let velocity_constraints = self.generate_velocity_constraints(&contacts);

        self.solve_velocity_constraints(contacts, &velocity_constraints);

        self.update_world();
    }

    fn apply_gravity(&mut self) {
        for body in self.rigid_bodies.values_mut() {
            body.apply_force_local(&(GRAVITY * body.mass), &Vector3::zeros());
        }
    }

    fn generate_contacts(&self) -> ArrayList<Contact> {
        let mut contacts = ArrayList::new();
        for ColliderPair(collider_id1, collider_id2) in self.broadphase.query_potential_collisions()
        {
            let collider1 =
                self.colliders.get(&collider_id1).expect("Erm... the collider aint here");
            let collider2 =
                self.colliders.get(&collider_id2).expect("Erm... the collider aint here");

            if collider1.body_id == collider2.body_id {
                continue;
            }

            /*
            let obb1 = self.obbs.get(&collider1.obb_id).expect("Collider has no OBB");
            let obb2 = self.obbs.get(&collider2.obb_id).expect("Collider has no OBB");

            if !obb1.collides(obb2) {
                continue;
            }
            */

            let (is_colliding, simplex) = gjk(collider1, collider2, &self.rigid_bodies);
            if is_colliding {
                let body1 = self
                    .rigid_bodies
                    .get(&collider1.body_id)
                    .expect("Collider doesn't belong to body");
                let body2 = self
                    .rigid_bodies
                    .get(&collider2.body_id)
                    .expect("Collider doesn't belong to body");

                if f32::is_nan(simplex[0].minkowski.x) {
                    continue;
                }

                match expanding_polytope(simplex, collider1, collider2, &self.rigid_bodies) {
                    Some((collision_normal, penetration_depth, collision_surface_points)) => {
                        let mean_collision_point =
                            (collision_surface_points.0 + collision_surface_points.1) / 2.;
                        contacts.push(Contact::build_contact(
                            body1,
                            body2,
                            mean_collision_point,
                            *collision_normal,
                            penetration_depth,
                        ));
                    }
                    None => continue,
                }
            }
        }
        return contacts;
    }

    fn generate_velocity_constraints(
        &self,
        contacts: &ArrayList<Contact>,
    ) -> ArrayList<VelocityConstraint> {
        let mut vel_constraints = ArrayList::with_capacity(3 * contacts.len());
        for contact in contacts {
            let (pen_jacobian, pen_bias) = contact.compute_normal_jacobian_bias(dt);
            vel_constraints.push(VelocityConstraint {
                body_id1: contact.body_id1,
                body_id2: contact.body_id2,

                jacobian: pen_jacobian,
                b: pen_bias,
            });

            let (fric_jacobian1, fric_jacobian2) = contact.compute_friction_jacobians();

            vel_constraints.push(VelocityConstraint {
                body_id1: contact.body_id1,
                body_id2: contact.body_id2,

                jacobian: fric_jacobian1,
                b: 0.,
            });

            vel_constraints.push(VelocityConstraint {
                body_id1: contact.body_id1,
                body_id2: contact.body_id2,

                jacobian: fric_jacobian2,
                b: 0.,
            });
        }

        return vel_constraints;
    }

    fn apply_impulse(&mut self, dv: Vector<f32, U12, ArrayStorage<f32, 12, 1>>, contact: &Contact) {
        let body1 =
            self.rigid_bodies.get_mut(&contact.body_id1).expect("Phantom body made collision");
        if !body1.is_static {
            body1.lin_velocity += Vector3::new(dv[0], dv[1], dv[2]);
            body1.ang_velocity += Vector3::new(dv[3], dv[4], dv[5]);
        }

        let body2 =
            self.rigid_bodies.get_mut(&contact.body_id2).expect("Phantom body made collision");
        if !body2.is_static {
            body2.lin_velocity += Vector3::new(dv[6], dv[7], dv[8]);
            body2.ang_velocity += Vector3::new(dv[9], dv[10], dv[11]);
        }
    }

    fn solve_velocity_constraints(
        &mut self,
        mut contacts: ArrayList<Contact>,
        constraints: &ArrayList<VelocityConstraint>,
    ) {
        for _ in 0..20 {
            for i in 0..contacts.len() {
                let contact = contacts.get_mut(i).expect("Array index out of bounds I guess");
                let normal_constraint =
                    constraints.get(3 * i).expect("Array index out of bounds I guess");
                let (delta_lambda, unscaled_velocity_response) =
                    normal_constraint.compute_constraint(&self.rigid_bodies);
                let tmp = contact.normal_impulse_sum;
                contact.normal_impulse_sum = f32::max(0., tmp + delta_lambda);
                let clamped_delta_lambda = contact.normal_impulse_sum - tmp;

                let delta_velocity = clamped_delta_lambda * unscaled_velocity_response;
                self.apply_impulse(delta_velocity, contact);

                let max_friction_magnitude = contact.coeff_friction * contact.normal_impulse_sum;

                let tangential_constraint1 =
                    constraints.get(3 * i + 1).expect("Array index out of bounds I guess");
                let (delta_lambda, unscaled_velocity_response) =
                    tangential_constraint1.compute_constraint(&self.rigid_bodies);
                let tmp = contact.tangent_impulse_sum1;
                contact.tangent_impulse_sum1 =
                    (tmp + delta_lambda).clamp(-max_friction_magnitude, max_friction_magnitude);
                let clamped_delta_lambda = contact.tangent_impulse_sum1 - tmp;

                let delta_velocity = clamped_delta_lambda * unscaled_velocity_response;
                self.apply_impulse(delta_velocity, contact);

                let tangential_constraint2 =
                    constraints.get(3 * i + 2).expect("Array index out of bounds I guess");
                let (delta_lambda, unscaled_velocity_response) =
                    tangential_constraint2.compute_constraint(&self.rigid_bodies);
                let tmp = contact.tangent_impulse_sum2;
                contact.tangent_impulse_sum2 =
                    (tmp + delta_lambda).clamp(-max_friction_magnitude, max_friction_magnitude);
                let clamped_delta_lambda = contact.tangent_impulse_sum2 - tmp;

                let delta_velocity = clamped_delta_lambda * unscaled_velocity_response;
                self.apply_impulse(delta_velocity, contact);
            }
        }
    }

    fn update_world(&mut self) {
        for rigid_body in self.rigid_bodies.values_mut() {
            rigid_body.update(dt);
        }

        for collider in self.colliders.values() {
            let obb = self.obbs.get_mut(&collider.obb_id).expect("Collider has no OBB");
            obb.update(&self.colliders, &self.rigid_bodies);
        }

        self.broadphase.update(&self.colliders, &self.obbs);
    }

    fn add_body(
        &mut self,
        position: Vector3<f32>,
        orientation: UnitQuaternion<f32>,
        initial_lin_velocity: Vector3<f32>,
        initial_ang_velocity: Vector3<f32>,
        collider_data: ArrayList<ColliderData>,
        is_static: bool,
    ) {
        // Collider datapoint: shape, mass, local center of mass, local orientation

        let body_id = self.next_body_id;
        let mut body = RigidBody::new_body(body_id, position, orientation, is_static);
        body.lin_velocity = initial_lin_velocity;
        body.ang_velocity = initial_ang_velocity;
        self.rigid_bodies.insert(body.id, body);
        self.next_body_id += 1;

        for data in collider_data {
            let body = self.rigid_bodies.get_mut(&body_id).expect("Body not added");

            let collider_id = self.next_collider_id;
            self.next_collider_id += 1;

            let obb_id = self.next_obb_id;
            self.next_obb_id += 1;

            let collider = Collider::new_collider(
                collider_id,
                body_id,
                obb_id,
                data.shape,
                data.mass,
                data.local_center_mass,
                data.local_orientation,
            );
            self.colliders.insert(collider.id, collider);
            body.add_collider(collider_id, &self.colliders);

            let obb = OBB::new(obb_id, collider_id, &self.colliders, &self.rigid_bodies);
            self.broadphase.add(AABB::new_from_obb(collider_id, &obb));
            self.obbs.insert(obb_id, obb);
        }
    }

    fn construct_world_box(&mut self, box_dimensions: [f32; 3]) {
        let bottom_wall_data = ColliderData {
            shape: ColliderShape::Box {
                x_len: (box_dimensions[0]),
                y_len: (box_dimensions[1]),
                z_len: (1.),
            },
            mass: f32::INFINITY,
            local_center_mass: Vector3::new(0., 0., -box_dimensions[2] / 2.),
            local_orientation: UnitQuaternion::identity(),
        };

        let top_wall_data = ColliderData {
            shape: ColliderShape::Box {
                x_len: (box_dimensions[0]),
                y_len: (box_dimensions[1]),
                z_len: (1.),
            },
            mass: f32::INFINITY,
            local_center_mass: Vector3::new(0., 0., box_dimensions[2] / 2.),
            local_orientation: UnitQuaternion::identity(),
        };

        let left_wall_data = ColliderData {
            shape: ColliderShape::Box {
                x_len: (box_dimensions[0]),
                y_len: (1.),
                z_len: (box_dimensions[2]),
            },
            mass: f32::INFINITY,
            local_center_mass: Vector3::new(0., -box_dimensions[1] / 2., 0.),
            local_orientation: UnitQuaternion::identity(),
        };

        let right_wall_data = ColliderData {
            shape: ColliderShape::Box {
                x_len: (box_dimensions[0]),
                y_len: (1.),
                z_len: (box_dimensions[2]),
            },
            mass: f32::INFINITY,
            local_center_mass: Vector3::new(0., box_dimensions[1] / 2., 0.),
            local_orientation: UnitQuaternion::identity(),
        };

        let front_wall_data = ColliderData {
            shape: ColliderShape::Box {
                x_len: (1.),
                y_len: (box_dimensions[1]),
                z_len: (box_dimensions[2]),
            },
            mass: f32::INFINITY,
            local_center_mass: Vector3::new(-box_dimensions[0] / 2., 0., 0.),
            local_orientation: UnitQuaternion::identity(),
        };

        let back_wall_data = ColliderData {
            shape: ColliderShape::Box {
                x_len: (1.),
                y_len: (box_dimensions[1]),
                z_len: (box_dimensions[2]),
            },
            mass: f32::INFINITY,
            local_center_mass: Vector3::new(box_dimensions[0] / 2., 0., 0.),
            local_orientation: UnitQuaternion::identity(),
        };

        self.add_body(
            Vector3::zeros(),
            UnitQuaternion::identity(),
            Vector3::zeros(),
            Vector3::zeros(),
            vec![
                bottom_wall_data,
                top_wall_data,
                left_wall_data,
                right_wall_data,
                front_wall_data,
                back_wall_data,
            ],
            true,
        );
    }

    pub fn initialize_world() -> World {
        let mut world = World {
            rigid_bodies: HashMap::new(),
            next_body_id: 0,
            colliders: HashMap::new(),
            next_collider_id: 0,
            obbs: HashMap::new(),
            next_obb_id: 0,
            broadphase: NSquared { aabb_list: ArrayList::new() },
        };

        world.construct_world_box([100., 100., 100.]);

        return world;
    }

    pub fn spawn_test_sphere(
        &mut self,
        radius: f32,
        mass: f32,
        position: Vector3<f32>,
        initial_lin_velocity: Vector3<f32>,
    ) {
        let sphere_data = ColliderData {
            shape: ColliderShape::Sphere { radius: radius },
            mass: mass,
            local_center_mass: Vector3::zeros(),
            local_orientation: UnitQuaternion::identity(),
        };

        self.add_body(
            position,
            UnitQuaternion::identity(),
            initial_lin_velocity,
            Vector3::zeros(),
            vec![sphere_data],
            false,
        );
    }

    pub fn spawn_test_box(
        &mut self,
        length: f32,
        width: f32,
        height: f32,
        mass: f32,
        position: Vector3<f32>,
        orientation: UnitQuaternion<f32>,
        initial_lin_velocity: Vector3<f32>,
        intitial_ang_velocity: Vector3<f32>,
    ) {
        let box_data = ColliderData {
            shape: ColliderShape::Box { x_len: length, y_len: width, z_len: height },
            mass: mass,
            local_center_mass: Vector3::zeros(),
            local_orientation: UnitQuaternion::identity(),
        };

        self.add_body(
            position,
            orientation,
            initial_lin_velocity,
            intitial_ang_velocity,
            vec![box_data],
            false,
        );
    }

    pub fn spawn_water_molecule(
        &mut self,
        o_radius: f32,
        o_mass: f32,
        bond_length: f32,
        position: Vector3<f32>,
        orientation: UnitQuaternion<f32>,
        initial_lin_velocity: Vector3<f32>,
        intitial_ang_velocity: Vector3<f32>,
    ) {
        let angle = 52.25_f32.to_radians();

        let h1_pos = Vector3::new(-bond_length * angle.cos(), -bond_length * angle.sin(), 0.);
        let h2_pos = Vector3::new(bond_length * angle.cos(), -bond_length * angle.sin(), 0.);

        let o_data = ColliderData {
            shape: ColliderShape::Sphere { radius: o_radius },
            mass: o_mass,
            local_center_mass: Vector3::zeros(),
            local_orientation: UnitQuaternion::identity(),
        };

        let h_data = |pos: Vector3<f32>| ColliderData {
            shape: ColliderShape::Sphere { radius: o_radius / 2. },
            mass: o_mass / 16.,
            local_center_mass: pos,
            local_orientation: UnitQuaternion::identity(),
        };

        let bond_data = |pos: Vector3<f32>| {
            let dir = pos.normalize();
            let rot = UnitQuaternion::rotation_between(&Vector3::x_axis(), &dir)
                .unwrap_or(UnitQuaternion::identity());
            ColliderData {
                shape: ColliderShape::Box {
                    x_len: bond_length,
                    y_len: o_radius / 6.,
                    z_len: o_radius / 6.,
                },
                mass: o_mass / 10000.,
                local_center_mass: pos / 2.,
                local_orientation: rot,
            }
        };

        self.add_body(
            position,
            orientation,
            initial_lin_velocity,
            intitial_ang_velocity,
            vec![o_data, h_data(h1_pos), h_data(h2_pos), bond_data(h1_pos), bond_data(h2_pos)],
            false,
        );
    }
}

struct ColliderData {
    shape: ColliderShape,
    mass: f32,
    local_center_mass: Vector3<f32>,
    local_orientation: UnitQuaternion<f32>,
}
