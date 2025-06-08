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
