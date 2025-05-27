use crate::{
    broadphase::ArrayList,
    rigid_body::{AABB, Collider, OBB, RigidBody},
    world::RigidBodyMap,
};
use core::{f32, panic};
use nalgebra::{Matrix3x2, SVD, UnitVector3, Vector3};

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
    vertices: [Vector3<f32>; 3],
    normal: UnitVector3<f32>,
    distance_from_origin: f32,
}

fn distance_origin_to_face(vertices: [Vector3<f32>; 3]) -> f32 {
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

fn build_face(vertices: [Vector3<f32>; 3]) -> Face {
    let [a, b, c] = vertices;

    let ab = b - a;
    let ac = c - a;
    let centroid = (a + b + c) / 3.;

    let mut norm = UnitVector3::new_normalize(ab.cross(&ac));
    if norm.dot(&centroid) < 0. {
        norm = -norm;
    }

    let distance = norm.dot(&a);

    return Face { vertices: vertices, normal: norm, distance_from_origin: distance };
}

fn build_initial_polytope(simplex: &ArrayList<Vector3<f32>>) -> ArrayList<Face> {
    let a = simplex[3];
    let b = simplex[2];
    let c = simplex[1];
    let d = simplex[0];

    return ArrayList::from([
        build_face([a, b, c]),
        build_face([a, b, d]),
        build_face([a, c, d]),
        build_face([b, c, d]),
    ]);
}

fn update_polytope(faces: ArrayList<Face>, new_point: Vector3<f32>) -> ArrayList<Face> {
    fn replace_face(face: &Face, new_point: Vector3<f32>) -> ArrayList<Face> {
        let mut new_faces = ArrayList::new();
        let [a, b, c] = face.vertices;
        new_faces.push(build_face([a, b, new_point]));
        new_faces.push(build_face([a, new_point, c]));
        new_faces.push(build_face([new_point, b, c]));

        return new_faces;
    }
    let mut new_faces: Vec<Face> = ArrayList::with_capacity(faces.len() * 3);

    for face in faces {
        let visible = (new_point - face.vertices[0]).dot(&face.normal) > 1e-6;
        if visible {
            new_faces.extend(replace_face(&face, new_point));
        } else {
            new_faces.push(face);
        }
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
    let mut faces = build_initial_polytope(&simplex);
    let mut closest_face = faces.get(0).expect("Faces list is empty");
    let mut min_distance = f32::INFINITY;
    loop {
        for face in faces {
            if face.distance_from_origin < min_distance {
                min_distance = face.distance_from_origin;
            }
        }
        let new_point = cso_support(collider1, collider2, &closest_face.normal, rigid_bodies).0;
        if new_point.norm() - min_distance < 1e-6 {
            break;
        }
        update_polytope(&mut faces, new_point);
    }
    return (closest_face.normal, closest_face.distance_from_origin);
}
