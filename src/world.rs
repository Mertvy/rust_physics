use nalgebra::{ArrayStorage, U12, Vector, Vector3};

use crate::broadphase::{ArrayList, BroadPhase, ColliderPair, NSquared};
use crate::collision::{Contact, expanding_polytope, gjk};
use crate::constraint::VelocityConstraint;
use crate::rigid_body::{self, Collider, RigidBody};
use std::collections::HashMap;

pub type RigidBodyMap = HashMap<usize, RigidBody>;
pub type ColliderMap = HashMap<usize, Collider>;

const dt: f32 = 1. / 60.;
const GRAVITY: Vector3<f32> = Vector3::new(0., 0., -9.8);

pub struct World {
    pub rigid_bodies: RigidBodyMap,
    pub colliders: ColliderMap,
    broadphase: NSquared,
}

impl World {
    fn step(&mut self) {
        self.apply_gravity();

        let contacts = self.generate_contacts();

        let velocity_constraints = self.generate_velocity_constraints(&contacts);

        self.solve_velocity_constraints(contacts, &velocity_constraints);

        self.update_world();
    }

    fn apply_gravity(&mut self) {
        for body in self.rigid_bodies.values_mut() {
            body.apply_force_local(&GRAVITY, &Vector3::zeros());
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

            let obb1 = &collider1.obb;
            let obb2 = &collider2.obb;

            if !obb1.collides(obb2) {
                continue;
            }

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

                let (collision_normal, penetration_depth, collision_surface_points) =
                    expanding_polytope(simplex, collider1, collider2, &self.rigid_bodies);
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
        body1.lin_velocity += Vector3::new(dv[0], dv[1], dv[2]);
        body1.ang_velocity += Vector3::new(dv[3], dv[4], dv[5]);

        let body2 =
            self.rigid_bodies.get_mut(&contact.body_id2).expect("Phantom body made collision");
        body2.lin_velocity += Vector3::new(dv[6], dv[7], dv[8]);
        body2.ang_velocity += Vector3::new(dv[9], dv[10], dv[11]);
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

        for collider in self.colliders.values_mut() {
            let rigid_body =
                self.rigid_bodies.get(&collider.body_id).expect("Collider does not belong to body");
            let mut obb = &collider.obb;
            obb.update(collider, rigid_body);
        }

        self.broadphase.update(&self.colliders);
    }
}
