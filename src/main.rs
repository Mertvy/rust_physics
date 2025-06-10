extern crate kiss3d;
extern crate nalgebra as na;

use kiss3d::light::Light;
use kiss3d::scene::SceneNode;
use kiss3d::window::Window;

use rayon;

use na::{Translation3, UnitQuaternion};
use std::collections::HashMap;

mod broadphase;
mod collision;
mod constraint;
mod rigid_body;
mod utils;
mod world;

use rigid_body::ColliderShape;
use world::World;

use crate::utils::rand_vector3;

fn main() {
    let mut window = Window::new("Physics Engine Renderer");
    window.set_light(Light::StickToCamera);

    let mut world = World::initialize_world([75., 75., 75.], 2000.);
    for _ in 0..150 {
        world.spawn_water_molecule(
            2.,
            1.,
            5.,
            rand_vector3(-35., 35.),
            UnitQuaternion::new(rand_vector3(-80., 80.)),
            rand_vector3(-30., 30.),
            rand_vector3(-10., 10.),
        );
        /*
        world.spawn_test_sphere(3., 1., rand_vector3(-50., 50.), rand_vector3(-5., 5.));
        world.spawn_test_box(
            10.,
            7.,
            5.,
            2.,
            rand_vector3(-50., 50.),
            UnitQuaternion::new(rand_vector3(-10., 10.)),
            rand_vector3(-5., 5.),
            rand_vector3(-0.5, 0.5),
        );
         */
    }

    // Track SceneNode per Collider ID
    let mut nodes: HashMap<usize, SceneNode> = HashMap::new();

    for (id, collider) in &world.colliders {
        let obb = world.obbs.get(&collider.obb_id).unwrap();
        let shape = &collider.shape;

        let mut node = match shape {
            ColliderShape::Box { x_len, y_len, z_len } => {
                let mut cube = window.add_cube(*x_len, *y_len, *z_len);
                cube.set_color(0.3, 0.6, 0.9);
                cube
            }
            ColliderShape::Sphere { radius } => {
                let mut sphere = window.add_sphere(*radius);
                sphere.set_color(0.9, 0.6, 0.3);
                sphere
            }
        };

        node.set_local_translation(Translation3::from(obb.global_pos));
        node.set_local_rotation(obb.global_orientation);

        nodes.insert(*id, node);
    }
    // Simulation/rendering loop
    while window.render() {
        world.step();

        for (id, collider) in &world.colliders {
            if let Some(node) = nodes.get_mut(id) {
                let obb = world.obbs.get(&collider.obb_id).unwrap();
                node.set_local_translation(Translation3::from(obb.global_pos));
                node.set_local_rotation(obb.global_orientation);
            }
        }
    }
}
