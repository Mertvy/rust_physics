use nalgebra::{ArrayStorage, Matrix3, RowVector, SMatrix, U12, Vector};

use crate::world::RigidBodyMap;

extern crate nalgebra as na;

pub struct VelocityConstraint {
    pub body_id1: usize,
    pub body_id2: usize,

    pub jacobian: RowVector<f32, U12, ArrayStorage<f32, 1, 12>>,
    pub b: f32,
}

fn block_diag_4x3x3(mats: [&Matrix3<f32>; 4]) -> SMatrix<f32, 12, 12> {
    let mut result = SMatrix::<f32, 12, 12>::zeros();

    result.fixed_view_mut::<3, 3>(0, 0).copy_from(mats[0]);
    result.fixed_view_mut::<3, 3>(3, 3).copy_from(mats[1]);
    result.fixed_view_mut::<3, 3>(6, 6).copy_from(mats[2]);
    result.fixed_view_mut::<3, 3>(9, 9).copy_from(mats[3]);

    return result;
}

impl VelocityConstraint {
    pub fn compute_constraint(
        &self,
        body_map: &RigidBodyMap,
    ) -> (f32, Vector<f32, U12, ArrayStorage<f32, 12, 1>>) {
        let body1 = body_map.get(&self.body_id1).expect("Collider does not exist");
        let body2 = body_map.get(&self.body_id2).expect("Collider does not exist");

        let velocities = na::stack![
            body1.lin_velocity;
            body1.ang_velocity;
            body2.lin_velocity;
            body2.ang_velocity
        ];

        let M_inv1 = Matrix3::identity() * body1.inv_mass;
        let M_inv2 = Matrix3::identity() * body2.inv_mass;
        let M_inv = block_diag_4x3x3([
            &M_inv1,
            &body1.global_inv_moment_inertia,
            &M_inv2,
            &body2.global_inv_moment_inertia,
        ]);

        let unscaled_velocity_reponse = M_inv * self.jacobian.transpose();
        let l = -((self.jacobian * velocities)[0] + self.b)
            / (self.jacobian * unscaled_velocity_reponse)[0];

        return (l, unscaled_velocity_reponse);
    }
}
