use crate::{
    broadphase::ArrayList,
    rigid_body::{AABB, Collider, OBB, RigidBody},
    world::RigidBodyMap,
};
use core::{f32, panic};
use nalgebra::{UnitVector3, Vector3};
use std::{
    collections::HashMap,
    f32::{INFINITY, consts::E},
    hash::Hash,
};

const THRESHOLD: f32 = 1e-6;

fn cso_support(
    collider1: &Collider,
    collider2: &Collider,
    dir: &Vector3<f32>,
    rigid_bodies: &RigidBodyMap,
) -> (Vector3<f32>, Vector3<f32>, Vector3<f32>) {
    let body1 = rigid_bodies.get(&collider1.body_id).expect("Collider has no body");
    let body2 = rigid_bodies.get(&collider2.body_id).expect("Collider has no body");

    let local_dir1 = body1.global_to_local_vec(&dir);
    let local_dir2 = body2.global_to_local_vec(&-dir);

    let support1 = collider1.support(&local_dir1);
    let support2 = collider2.support(&local_dir2);

    let support1 = body1.local_to_global_pos(&support1);
    let support2 = body2.local_to_global_pos(&support2);

    return (support1 - support2, support1, support2);
}

enum SimplexUpdate {
    ContainsOrigin,
    Continue { simplex: ArrayList<Vector3<f32>>, direction: UnitVector3<f32> },
}

fn update_simplex(simplex: &mut ArrayList<Vector3<f32>>) -> SimplexUpdate {
    fn update_line(simplex: &mut ArrayList<Vector3<f32>>) -> SimplexUpdate {
        let a = simplex[1];
        let b = simplex[0];
        let ab = b - a;
        let ao = -a;

        if ab.dot(&ao) >= 0. {
            let ab_perp = ab.cross(&ao).cross(&ab);
            return SimplexUpdate::Continue {
                simplex: vec![b, a],
                direction: UnitVector3::new_normalize(ab_perp),
            };
        } else {
            return SimplexUpdate::Continue {
                simplex: vec![a],
                direction: UnitVector3::new_normalize(ao),
            };
        }
    }

    fn update_triangle(simplex: &mut ArrayList<Vector3<f32>>) -> SimplexUpdate {
        let a = simplex[2];
        let b = simplex[1];
        let c = simplex[0];

        let ab = b - a;
        let ac = c - a;
        let ao = -a;

        let abc = ab.cross(&ac);

        let ab_perp = abc.cross(&ab);
        if ab_perp.dot(&ao) > 0.0 {
            return SimplexUpdate::Continue {
                simplex: vec![b, a],
                direction: UnitVector3::new_normalize(ab.cross(&ao).cross(&ab)),
            };
        }

        let ac_perp = ac.cross(&abc);
        if ac_perp.dot(&ao) > 0.0 {
            return SimplexUpdate::Continue {
                simplex: vec![c, a],
                direction: UnitVector3::new_normalize(ac.cross(&ao).cross(&ac)),
            };
        }

        // Origin is in triangle region
        let dir = if abc.dot(&ao) > 0.0 { abc } else { -abc };
        return SimplexUpdate::Continue {
            simplex: vec![c, b, a],
            direction: UnitVector3::new_normalize(dir),
        };
    }

    fn update_tetrahedron(simplex: &mut ArrayList<Vector3<f32>>) -> SimplexUpdate {
        let a = simplex[3];
        let b = simplex[2];
        let c = simplex[1];
        let d = simplex[0];
        let ao = -a;

        let ab = b - a;
        let ac = c - a;
        let ad = d - a;

        let mut abc = ab.cross(&ac);
        if abc.dot(&ad) > 0.0 {
            abc = -abc;
        }

        let mut acd = ac.cross(&ad);
        if acd.dot(&ab) > 0.0 {
            acd = -acd;
        }

        let mut adb = ad.cross(&ab);
        if adb.dot(&ac) > 0.0 {
            adb = -adb;
        }

        let inside_abc = abc.dot(&ao) > THRESHOLD;
        let inside_acd = acd.dot(&ao) > THRESHOLD;
        let inside_adb = adb.dot(&ao) > THRESHOLD;

        if inside_abc {
            return SimplexUpdate::Continue {
                simplex: vec![c, b, a],
                direction: UnitVector3::new_normalize(abc),
            };
        }

        if inside_acd {
            return SimplexUpdate::Continue {
                simplex: vec![d, c, a],
                direction: UnitVector3::new_normalize(acd),
            };
        }

        if inside_adb {
            return SimplexUpdate::Continue {
                simplex: vec![b, d, a],
                direction: UnitVector3::new_normalize(adb),
            };
        }

        return SimplexUpdate::ContainsOrigin;
    }

    match simplex.len() {
        2 => update_line(simplex),
        3 => update_triangle(simplex),
        4 => update_tetrahedron(simplex),
        _ => panic!("Invalid simplex dimension"),
    }
}

fn gjk(
    collider1: &Collider,
    collider2: &Collider,
    rigid_bodies: &RigidBodyMap,
) -> (bool, ArrayList<Vector3<f32>>) {
    // https://en.wikipedia.org/wiki/Gilbert%E2%80%93Johnson%E2%80%93Keerthi_distance_algorithm

    let mut simplex = ArrayList::new();
    let mut dir = UnitVector3::new_normalize(Vector3::new(1., 0., 0.));
    simplex.push(cso_support(collider1, collider2, &dir, rigid_bodies).0);
    if simplex[0].dot(&dir) <= 0. {
        return (false, simplex);
    }
    dir = UnitVector3::new_normalize(-simplex[0]);
    loop {
        let new_point = cso_support(collider1, collider2, &dir, rigid_bodies).0;
        for point in &simplex {
            if (point - new_point).norm_squared() < THRESHOLD {
                return (false, simplex);
            }
        }
        simplex.push(new_point);
        if new_point.dot(&dir) <= 0. {
            return (false, simplex);
        }

        let update = update_simplex(&mut simplex);
        match update {
            SimplexUpdate::ContainsOrigin => return (true, simplex),
            SimplexUpdate::Continue { simplex: new_simplex, direction: new_dir } => {
                dir = new_dir;
                simplex.clear();
                simplex.extend(new_simplex);
            }
        }
    }
}

fn distance_origin_to_triangle(vertices: [Vector3<f32>; 3]) -> f32 {
    let [a, b, c] = vertices;

    let ab = b - a;
    let ac = c - a;
    let normal = ab.cross(&ac).normalize();

    let dist = normal.dot(&a);
    let proj = -dist * normal;

    let ap = proj - a;
    let bp = proj - b;
    let cp = proj - c;

    let edge0 = b - a;
    let edge1 = c - b;
    let edge2 = a - c;

    let c0 = edge0.cross(&ap);
    let c1 = edge1.cross(&bp);
    let c2 = edge2.cross(&cp);

    let inside = c0.dot(&normal) >= 0.0 && c1.dot(&normal) >= 0.0 && c2.dot(&normal) >= 0.0;

    if inside {
        dist.abs()
    } else {
        fn point_segment_distance(p: Vector3<f32>, a: Vector3<f32>, b: Vector3<f32>) -> f32 {
            let ab = b - a;
            let t = (p - a).dot(&ab) / ab.norm_squared();
            let t_clamped = t.clamp(0.0, 1.0);
            let projection = a + t_clamped * ab;
            (p - projection).norm()
        }

        let d0 = point_segment_distance(Vector3::zeros(), a, b);
        let d1 = point_segment_distance(Vector3::zeros(), b, c);
        let d2 = point_segment_distance(Vector3::zeros(), c, a);
        d0.min(d1).min(d2)
    }
}

fn construct_orthonormal_basis(v1: Vector3<f32>) -> [UnitVector3<f32>; 3] {
    let v1 = UnitVector3::new_normalize(v1);
    let v2 = if v1.x >= 0.57735 {
        UnitVector3::new_normalize(Vector3::new(v1.y, -v1.x, 0.))
    } else {
        UnitVector3::new_normalize(Vector3::new(0., v1.z, -v1.y))
    };
    let v3 = UnitVector3::new_normalize(v1.cross(&v2));

    return [v1, v2, v3];
}

fn expand_simplex(
    simplex: &mut ArrayList<Vector3<f32>>,
    collider1: &Collider,
    collider2: &Collider,
    rigid_bodies: &RigidBodyMap,
) {
    if simplex.len() == 2 {
        let basis = construct_orthonormal_basis(simplex[0] - simplex[1]);
        simplex.push(cso_support(collider1, collider2, &basis[1], rigid_bodies).0);
        simplex.push(cso_support(collider1, collider2, &basis[2], rigid_bodies).0);
    } else if simplex.len() == 3 {
        let a = simplex[2];
        let b = simplex[1];
        let c = simplex[0];
        let dir = UnitVector3::new_normalize((b - a).cross(&(c - a)));
        simplex.push(cso_support(collider1, collider2, &dir, rigid_bodies).0);
    }
}

struct Face {
    vertex_ids: [usize; 3],
    normal: UnitVector3<f32>,
    distance_from_origin: f32,
}

impl Face {
    fn build_face(vertex_ids: [usize; 3], vertex_map: &HashMap<usize, Vector3<f32>>) -> Face {
        let [a_id, b_id, c_id] = vertex_ids;
        let a = vertex_map.get(&a_id).expect("Vertex does not exist");
        let b = vertex_map.get(&b_id).expect("Vertex does not exist");
        let c = vertex_map.get(&c_id).expect("Vertex does not exist");

        let ab = b - a;
        let ac = c - a;
        let centroid = (a + b + c) / 3.;

        let mut norm = UnitVector3::new_normalize(ab.cross(&ac));
        if norm.dot(&centroid) < 0. {
            norm = -norm;
        }

        let distance = norm.dot(&a);

        return Face { vertex_ids: vertex_ids, normal: norm, distance_from_origin: distance };
    }
}

fn build_initial_polytope(
    vertex_ids: [usize; 4],
    vertex_map: &HashMap<usize, Vector3<f32>>,
) -> ArrayList<Face> {
    let [d_id, c_id, b_id, a_id] = vertex_ids;

    return ArrayList::from([
        Face::build_face([a_id, b_id, c_id], vertex_map),
        Face::build_face([a_id, b_id, d_id], vertex_map),
        Face::build_face([a_id, c_id, d_id], vertex_map),
        Face::build_face([b_id, c_id, d_id], vertex_map),
    ]);
}

#[derive(Hash, PartialEq, std::cmp::Eq)]
struct Edge(usize, usize);

impl Edge {
    fn new(v1: usize, v2: usize) -> Edge {
        return Edge(usize::min(v1, v2), usize::max(v1, v2));
    }
}

fn update_polytope(
    faces: ArrayList<Face>,
    new_vertex_id: usize,
    vertex_map: &HashMap<usize, Vector3<f32>>,
) -> ArrayList<Face> {
    let mut new_faces = ArrayList::new();

    let new_vertex = vertex_map.get(&new_vertex_id).expect("New point does not exist");
    let mut visible_edges_count: HashMap<Edge, usize> = HashMap::new();

    for face in faces {
        let face_corner = vertex_map.get(&face.vertex_ids[0]).expect("Face vertex missing");
        let visible = (new_vertex - face_corner).dot(&face.normal) > 1e-6;

        if visible {
            let edges = [
                Edge::new(face.vertex_ids[0], face.vertex_ids[1]),
                Edge::new(face.vertex_ids[0], face.vertex_ids[2]),
                Edge::new(face.vertex_ids[1], face.vertex_ids[2]),
            ];

            for edge in edges {
                *visible_edges_count.entry(edge).or_insert(0) += 1;
            }
        } else {
            new_faces.push(face);
        }
    }

    let mut boundary_edges = ArrayList::new();
    for (edge, count) in visible_edges_count {
        if count == 1 {
            boundary_edges.push(edge);
        }
    }

    for edge in boundary_edges {
        let new_face_vertex_ids = [edge.0, edge.1, new_vertex_id];
        new_faces.push(Face::build_face(new_face_vertex_ids, vertex_map));
    }

    return new_faces;
}

fn expanding_polytope(
    mut simplex: ArrayList<Vector3<f32>>,
    collider1: &Collider,
    collider2: &Collider,
    rigid_bodies: &RigidBodyMap,
) -> (UnitVector3<f32>, f32) {
    expand_simplex(&mut simplex, collider1, collider2, rigid_bodies);

    let mut vertex_map: HashMap<usize, Vector3<f32>> = HashMap::new();
    for i in 0..4 {
        vertex_map.insert(i, simplex[i]);
    }
    let mut new_vertex_id = 3;

    let mut faces = build_initial_polytope([0, 1, 2, 3], &vertex_map);
    let mut closest_face;

    loop {
        closest_face = faces
            .iter()
            .min_by(|a, b| a.distance_from_origin.partial_cmp(&b.distance_from_origin).unwrap())
            .unwrap();

        let new_vertex = cso_support(collider1, collider2, &closest_face.normal, rigid_bodies).0;
        new_vertex_id += 1;
        vertex_map.insert(new_vertex_id, new_vertex);

        let projected_distance = new_vertex.dot(&closest_face.normal);
        if projected_distance - closest_face.distance_from_origin < 1e-6 {
            break;
        }
        faces = update_polytope(faces, new_vertex_id, &vertex_map);
    }
    return (closest_face.normal, closest_face.distance_from_origin);
}
