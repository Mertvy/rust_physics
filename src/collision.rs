use crate::{
    broadphase::ArrayList,
    rigid_body::{Collider, RigidBody},
    utils::safe_normalize,
    world::RigidBodyMap,
};
use core::{f32, panic};
use nalgebra::{ArrayStorage, RowVector, U12, UnitVector3, Vector3};
use std::{collections::HashMap, hash::Hash};

const THRESHOLD: f32 = 1e-6;

#[derive(Clone, Copy)]
pub struct SupportPoint {
    pub minkowski: Vector3<f32>,
    support1: Vector3<f32>,
    support2: Vector3<f32>,
}

fn cso_support(
    collider1: &Collider,
    collider2: &Collider,
    dir: &Vector3<f32>,
    rigid_bodies: &RigidBodyMap,
) -> SupportPoint {
    let body1 = rigid_bodies.get(&collider1.body_id).expect("Collider has no body");
    let body2 = rigid_bodies.get(&collider2.body_id).expect("Collider has no body");

    let local_dir1 = body1.global_to_local_vec(&dir);
    let local_dir2 = body2.global_to_local_vec(&-dir);

    let support1 = collider1.support(&local_dir1);
    let support2 = collider2.support(&local_dir2);

    let support1 = body1.local_to_global_pos(&support1);
    let support2 = body2.local_to_global_pos(&support2);

    return SupportPoint { minkowski: support1 - support2, support1: support1, support2: support2 };
}

enum SimplexUpdate {
    ContainsOrigin,
    Continue { simplex: ArrayList<SupportPoint>, direction: UnitVector3<f32> },
}

fn update_simplex(simplex: &mut ArrayList<SupportPoint>) -> SimplexUpdate {
    fn update_line(simplex: &mut ArrayList<SupportPoint>) -> SimplexUpdate {
        let a = simplex[1];
        let b = simplex[0];
        let ab = b.minkowski - a.minkowski;
        let ao = -a.minkowski;

        if ab.dot(&ao) < 0. {
            return SimplexUpdate::Continue { simplex: vec![a], direction: safe_normalize(ao) };
        }

        let ab_perp = ab.cross(&ao).cross(&ab);
        return SimplexUpdate::Continue { simplex: vec![b, a], direction: safe_normalize(ab_perp) };
    }

    fn update_triangle(simplex: &mut ArrayList<SupportPoint>) -> SimplexUpdate {
        let a = simplex[2];
        let b = simplex[1];
        let c = simplex[0];

        let ab = b.minkowski - a.minkowski;
        let ac = c.minkowski - a.minkowski;
        let ao = -a.minkowski;

        let abc = ab.cross(&ac);

        let mut ab_perp = abc.cross(&ab);
        if ab_perp.dot(&ac) > 0. {
            ab_perp = -ab_perp;
        }
        if ab_perp.dot(&ao) > 0. {
            return update_line(&mut vec![b, a]);
        }

        let mut ac_perp = ac.cross(&abc);
        if ac_perp.dot(&ab) > 0. {
            ac_perp = -ac_perp;
        }
        if ac_perp.dot(&ao) > 0. {
            return update_line(&mut vec![c, a]);
        }

        // Origin is in triangle region
        let dir = if abc.dot(&ao) > 0. { abc } else { -abc };
        return SimplexUpdate::Continue { simplex: vec![c, b, a], direction: safe_normalize(dir) };
    }

    fn update_tetrahedron(simplex: &mut ArrayList<SupportPoint>) -> SimplexUpdate {
        let a = simplex[3];
        let b = simplex[2];
        let c = simplex[1];
        let d = simplex[0];
        let ao = -a.minkowski;

        let ab = b.minkowski - a.minkowski;
        let ac = c.minkowski - a.minkowski;
        let ad = d.minkowski - a.minkowski;

        let mut abc = ab.cross(&ac);
        if abc.dot(&ad) > 0. {
            abc = -abc;
        }

        let mut acd = ac.cross(&ad);
        if acd.dot(&ab) > 0. {
            acd = -acd;
        }

        let mut adb = ad.cross(&ab);
        if adb.dot(&ac) > 0. {
            adb = -adb;
        }

        let inside_abc = abc.dot(&ao) > THRESHOLD;
        let inside_acd = acd.dot(&ao) > THRESHOLD;
        let inside_adb = adb.dot(&ao) > THRESHOLD;

        if inside_abc {
            return SimplexUpdate::Continue {
                simplex: vec![c, b, a],
                direction: safe_normalize(abc),
            };
        }

        if inside_acd {
            return SimplexUpdate::Continue {
                simplex: vec![d, c, a],
                direction: safe_normalize(acd),
            };
        }

        if inside_adb {
            return SimplexUpdate::Continue {
                simplex: vec![b, d, a],
                direction: safe_normalize(adb),
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

pub fn gjk(
    collider1: &Collider,
    collider2: &Collider,
    rigid_bodies: &RigidBodyMap,
) -> (bool, ArrayList<SupportPoint>) {
    // https://en.wikipedia.org/wiki/Gilbert%E2%80%93Johnson%E2%80%93Keerthi_distance_algorithm

    let mut simplex = ArrayList::new();
    let mut dir = safe_normalize(Vector3::new(1., 0., 0.));
    simplex.push(cso_support(collider1, collider2, &dir, rigid_bodies));
    if simplex[0].minkowski.dot(&dir) <= 0. {
        return (false, simplex);
    }
    dir = safe_normalize(-simplex[0].minkowski);

    loop {
        let new_support = cso_support(collider1, collider2, &dir, rigid_bodies);
        for point in &simplex {
            if (point.minkowski - new_support.minkowski).norm_squared() < THRESHOLD {
                return (false, simplex);
            }
        }
        simplex.push(new_support);
        if new_support.minkowski.dot(&dir) <= 0. {
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

    let inside = c0.dot(&normal) >= 0. && c1.dot(&normal) >= 0. && c2.dot(&normal) >= 0.;

    if inside {
        dist.abs()
    } else {
        fn point_segment_distance(p: Vector3<f32>, a: Vector3<f32>, b: Vector3<f32>) -> f32 {
            let ab = b - a;
            let t = (p - a).dot(&ab) / ab.norm_squared();
            let t_clamped = t.clamp(0., 1.0);
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
    let v1 = safe_normalize(v1);
    let v2 = if v1.x >= 0.57735 {
        safe_normalize(Vector3::new(v1.y, -v1.x, 0.))
    } else {
        safe_normalize(Vector3::new(0., v1.z, -v1.y))
    };
    let v3 = safe_normalize(v1.cross(&v2));

    return [v1, v2, v3];
}

fn expand_simplex(
    simplex: &mut ArrayList<SupportPoint>,
    collider1: &Collider,
    collider2: &Collider,
    rigid_bodies: &RigidBodyMap,
) {
    if simplex.len() == 2 {
        let basis = construct_orthonormal_basis(simplex[0].minkowski - simplex[1].minkowski);
        simplex.push(cso_support(collider1, collider2, &basis[1], rigid_bodies));
        simplex.push(cso_support(collider1, collider2, &basis[2], rigid_bodies));
    } else if simplex.len() == 3 {
        let a = simplex[2];
        let b = simplex[1];
        let c = simplex[0];
        let dir = safe_normalize((b.minkowski - a.minkowski).cross(&(c.minkowski - a.minkowski)));
        simplex.push(cso_support(collider1, collider2, &dir, rigid_bodies));
    }
}

#[derive(Debug)]
struct Face {
    vertex_ids: [usize; 3],
    normal: UnitVector3<f32>,
    distance_from_origin: f32,
}

impl Face {
    fn build_face(vertex_ids: [usize; 3], vertex_map: &HashMap<usize, SupportPoint>) -> Face {
        let [a_id, b_id, c_id] = vertex_ids;
        let a = vertex_map.get(&a_id).expect("Vertex does not exist");
        let b = vertex_map.get(&b_id).expect("Vertex does not exist");
        let c = vertex_map.get(&c_id).expect("Vertex does not exist");

        let ab = b.minkowski - a.minkowski;
        let ac = c.minkowski - a.minkowski;
        let centroid = (a.minkowski + b.minkowski + c.minkowski) / 3.;

        let mut norm = safe_normalize(ab.cross(&ac));
        if norm.dot(&centroid) < 0. {
            norm = -norm;
        }

        let distance = norm.dot(&a.minkowski);

        /*
        eprintln!(
            "Building face from: {:?}, face normal = {:?}, distance from origin = {:?}",
            vertex_ids, norm, distance
        );
         */
        return Face { vertex_ids: vertex_ids, normal: norm, distance_from_origin: distance };
    }
}

fn build_initial_polytope(
    vertex_ids: [usize; 4],
    vertex_map: &HashMap<usize, SupportPoint>,
) -> ArrayList<Face> {
    let [d_id, c_id, b_id, a_id] = vertex_ids;

    return ArrayList::from([
        Face::build_face([a_id, b_id, c_id], vertex_map),
        Face::build_face([a_id, b_id, d_id], vertex_map),
        Face::build_face([a_id, c_id, d_id], vertex_map),
        Face::build_face([b_id, c_id, d_id], vertex_map),
    ]);
}

#[derive(Debug, Hash, PartialEq, Eq, Clone, Copy)]
struct Edge(usize, usize);

impl Edge {
    fn new(a: usize, b: usize) -> Self {
        if a < b { Edge(a, b) } else { Edge(b, a) }
    }
}

struct Polytope {
    faces: Vec<Face>,
    active_faces: Vec<usize>, // IDs into faces
}

impl Polytope {
    fn build_initial(vertex_ids: [usize; 4], vertex_map: &HashMap<usize, SupportPoint>) -> Self {
        let [d, c, b, a] = vertex_ids;

        let mut faces = vec![
            Face::build_face([a, b, c], vertex_map),
            Face::build_face([a, b, d], vertex_map),
            Face::build_face([a, c, d], vertex_map),
            Face::build_face([b, c, d], vertex_map),
        ];

        let active_faces = vec![0, 1, 2, 3];
        Self { faces, active_faces }
    }

    fn get_closest_face(&self) -> usize {
        *self
            .active_faces
            .iter()
            .min_by(|&&a, &&b| {
                self.faces[a]
                    .distance_from_origin
                    .partial_cmp(&self.faces[b].distance_from_origin)
                    .unwrap()
            })
            .expect("Polytope has no faces")
    }

    fn expand(&mut self, new_vertex_id: usize, vertex_map: &HashMap<usize, SupportPoint>) -> bool {
        let new_vertex = &vertex_map[&new_vertex_id];

        let mut visible_faces = vec![];
        let mut edge_counts: HashMap<Edge, usize> = HashMap::new();
        let mut any_faces_visible = false;

        self.active_faces.retain(|&face_id| {
            let face = &self.faces[face_id];
            let face_corner = &vertex_map[&face.vertex_ids[0]];
            let visible = (new_vertex.minkowski - face_corner.minkowski).dot(&face.normal) > 1e-6;

            if visible {
                any_faces_visible = true;
                visible_faces.push(face_id);

                let [a, b, c] = face.vertex_ids;
                for edge in [Edge::new(a, b), Edge::new(a, c), Edge::new(b, c)] {
                    *edge_counts.entry(edge).or_insert(0) += 1;
                }

                false // remove this face
            } else {
                true // keep this face
            }
        });

        if !any_faces_visible {
            return true;
        }

        let boundary_edges: Vec<_> = edge_counts
            .into_iter()
            .filter_map(|(e, count)| if count == 1 { Some(e) } else { None })
            .collect();

        if boundary_edges.is_empty() {
            return false;
        }

        for edge in boundary_edges {
            let ids = [edge.0, edge.1, new_vertex_id];
            let face = Face::build_face(ids, vertex_map);
            let face_id = self.faces.len();
            self.faces.push(face);
            self.active_faces.push(face_id);
        }
        return true;
    }
}

fn update_polytope(
    faces: ArrayList<Face>,
    new_vertex_id: usize,
    vertex_map: &HashMap<usize, SupportPoint>,
) -> ArrayList<Face> {
    let mut new_faces = ArrayList::new();

    let new_vertex = vertex_map.get(&new_vertex_id).expect("New point does not exist");
    let mut visible_edges_count: HashMap<Edge, usize> = HashMap::new();

    for face in faces {
        let face_corner = vertex_map.get(&face.vertex_ids[0]).expect("Face vertex missing");
        let visible = (new_vertex.minkowski - face_corner.minkowski).dot(&face.normal) > 1e-6;

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

    if boundary_edges.is_empty() {
        panic!("EPA failed: No boundary edges found after expanding polytope");
    }

    for edge in boundary_edges {
        let new_face_vertex_ids = [edge.0, edge.1, new_vertex_id];
        new_faces.push(Face::build_face(new_face_vertex_ids, vertex_map));
    }

    return new_faces;
}

fn compute_suface_collision_points(
    face: &Face,
    vertex_map: &HashMap<usize, SupportPoint>,
) -> (Vector3<f32>, Vector3<f32>) {
    let a = vertex_map.get(&face.vertex_ids[0]).expect("Face doesnt have vertices");
    let b = vertex_map.get(&face.vertex_ids[1]).expect("Face doesnt have vertices");
    let c = vertex_map.get(&face.vertex_ids[2]).expect("Face doesnt have vertices");

    let p = face.normal.into_inner() * face.distance_from_origin;

    let v0 = b.minkowski - a.minkowski;
    let v1 = c.minkowski - a.minkowski;
    let v2 = p - a.minkowski;

    let d00 = v0.dot(&v0);
    let d01 = v0.dot(&v1);
    let d11 = v1.dot(&v1);
    let d20 = v2.dot(&v0);
    let d21 = v2.dot(&v1);

    let denom = d00 * d11 - d01 * d01;
    let v = (d11 * d20 - d01 * d21) / denom;
    let w = (d00 * d21 - d01 * d20) / denom;
    let u = 1.0 - v - w;

    if denom.abs() < 1e-8 {
        // Degenerate triangle, return something safe or fallback
        let avg1 = (a.support1 + b.support1 + c.support1) / 3.0;
        let avg2 = (a.support2 + b.support2 + c.support2) / 3.0;
        return (avg1, avg2);
    }

    return (
        u * a.support1 + v * b.support1 + w * c.support1,
        u * a.support2 + v * b.support2 + w * c.support2,
    );
}

pub fn expanding_polytope(
    mut simplex: ArrayList<SupportPoint>,
    collider1: &Collider,
    collider2: &Collider,
    rigid_bodies: &RigidBodyMap,
) -> Option<(UnitVector3<f32>, f32, (Vector3<f32>, Vector3<f32>))> {
    expand_simplex(&mut simplex, collider1, collider2, rigid_bodies);

    let mut vertex_map = HashMap::new();
    for (i, p) in simplex.into_iter().enumerate() {
        vertex_map.insert(i, p);
    }
    let mut new_vertex_id = vertex_map.len();

    let mut polytope = Polytope::build_initial([0, 1, 2, 3], &vertex_map);

    let mut i = 0;
    const MAX_ITER: usize = 64;
    loop {
        let best_face_id = polytope.get_closest_face();
        let best_face = &polytope.faces[best_face_id];

        let new_point = cso_support(collider1, collider2, &best_face.normal, rigid_bodies);
        vertex_map.insert(new_vertex_id, new_point);

        let projected_dist = new_point.minkowski.dot(&best_face.normal);
        if (projected_dist - best_face.distance_from_origin).abs() < 1e-6 {
            let collision_points = compute_suface_collision_points(best_face, &vertex_map);
            return Some((best_face.normal, best_face.distance_from_origin, collision_points));
        }

        if i >= MAX_ITER {
            let collision_points = compute_suface_collision_points(best_face, &vertex_map);
            return Some((best_face.normal, best_face.distance_from_origin, collision_points));
        }
        i += 1;

        if !polytope.expand(new_vertex_id, &vertex_map) {
            return None;
        }
        new_vertex_id += 1;
    }
}
pub struct Contact {
    pub body_id1: usize,
    pub body_id2: usize,

    com1: Vector3<f32>,
    com2: Vector3<f32>,

    lin_velocity1: Vector3<f32>,
    ang_velocity1: Vector3<f32>,

    lin_velocity2: Vector3<f32>,
    ang_velocity2: Vector3<f32>,

    mean_collision_point: Vector3<f32>,

    normal: Vector3<f32>, // from collider 1 to 2
    tangent1: Vector3<f32>,
    tangent2: Vector3<f32>,

    penetration_depth: f32,

    pub normal_impulse_sum: f32,
    pub tangent_impulse_sum1: f32,
    pub tangent_impulse_sum2: f32,

    pub coeff_friction: f32,
}

fn build_jacobian(vectors: [Vector3<f32>; 4]) -> RowVector<f32, U12, ArrayStorage<f32, 1, 12>> {
    let [a, b, c, d] = vectors;
    return RowVector::<f32, U12, ArrayStorage<f32, 1, 12>>::from_row_slice(&[
        a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z, d.x, d.y, d.z,
    ]);
}

impl Contact {
    pub fn build_contact(
        body1: &RigidBody,
        body2: &RigidBody,
        mean_collision_point: Vector3<f32>,
        normal: Vector3<f32>,
        penetration_depth: f32,
    ) -> Contact {
        let [normal, tangent1, tangent2] = construct_orthonormal_basis(normal);
        return Contact {
            body_id1: body1.id,
            body_id2: body2.id,

            com1: body1.global_center_mass,
            com2: body2.global_center_mass,

            lin_velocity1: body1.lin_velocity,
            ang_velocity1: body1.ang_velocity,

            lin_velocity2: body2.lin_velocity,
            ang_velocity2: body2.ang_velocity,

            mean_collision_point: mean_collision_point,

            normal: *normal,
            tangent1: *tangent1,
            tangent2: *tangent2,

            penetration_depth: penetration_depth,

            normal_impulse_sum: 0.,
            tangent_impulse_sum1: 0.,
            tangent_impulse_sum2: 0.,

            coeff_friction: 0.2,
        };
    }

    pub fn compute_normal_jacobian_bias(
        &self,
        dt: f32,
    ) -> (RowVector<f32, U12, ArrayStorage<f32, 1, 12>>, f32) {
        let r1 = self.mean_collision_point - self.com1;
        let r2 = self.mean_collision_point - self.com2;
        let jacobian = build_jacobian([
            -self.normal,
            -r1.cross(&self.normal),
            self.normal,
            r2.cross(&self.normal),
        ]);
        let beta = 0.1; // Tunable
        let baumgarte_stabilization = -beta * self.penetration_depth / dt;

        let rel_velocity = (-self.lin_velocity1 - self.ang_velocity1.cross(&r1)
            + self.lin_velocity2
            + self.ang_velocity2.cross(&r2))
        .dot(&self.normal);

        let restitution_coeff = 0.5; // Tunable 
        let restitution_bias = if rel_velocity < -1e-2 && self.penetration_depth > 1e-3 {
            restitution_coeff
                * (-self.lin_velocity1 - self.ang_velocity1.cross(&r1)
                    + self.lin_velocity2
                    + self.ang_velocity2.cross(&r2))
                .dot(&self.normal)
        } else {
            0.
        };

        let b = baumgarte_stabilization + restitution_bias;
        return (jacobian, b);
    }
    pub fn compute_friction_jacobians(
        &self,
    ) -> (
        RowVector<f32, U12, ArrayStorage<f32, 1, 12>>,
        RowVector<f32, U12, ArrayStorage<f32, 1, 12>>,
    ) {
        let r1 = self.mean_collision_point - self.com1;
        let r2 = self.mean_collision_point - self.com2;
        let jacobian1 = build_jacobian([
            -self.tangent1,
            -r1.cross(&self.tangent1),
            self.tangent1,
            r2.cross(&self.tangent1),
        ]);
        let jacobian2 = build_jacobian([
            -self.tangent2,
            -r1.cross(&self.tangent2),
            self.tangent2,
            r2.cross(&self.tangent2),
        ]);
        return (jacobian1, jacobian2);
    }
}
