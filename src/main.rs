extern crate kiss3d;
extern crate nalgebra as na;

use kiss3d::light::Light;
use kiss3d::scene::SceneNode;
use kiss3d::window::Window;

use na::{Translation3, UnitQuaternion, Vector3};
use std::collections::HashMap;

mod broadphase;
mod collision;
mod constraint;
mod rigid_body;
mod utils;
mod world;

use rigid_body::{ColliderShape, OBB};
use world::World;

use crate::utils::rand_vector3;

fn main() {
    let mut window = Window::new("Physics Engine Renderer");
    window.set_light(Light::StickToCamera);

    let mut world = World::initialize_world();
    for _ in 0..50 {
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
    }

    println!("Number of colliders: {}", world.colliders.len());
    for (id, collider) in &world.colliders {
        println!("Collider {id}: {:?}", collider.shape);
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

    println!("Total SceneNodes: {}", nodes.len());
    for (id, node) in &nodes {
        println!("Node {id} exists");
    }

    println!("Collider 9 position: {:?}", world.obbs.get(&9).unwrap().global_pos);
    println!("Collider 9 orientation: {:?}", world.obbs.get(&9).unwrap().global_orientation);
    println!("Collider 10 position: {:?}", world.obbs.get(&10).unwrap().global_pos);
    println!("Collider 10 orientation: {:?}", world.obbs.get(&10).unwrap().global_orientation);

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
