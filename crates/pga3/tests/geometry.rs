//! Euclidean geometry oracles for the PGA embedding.

use pga3::{
    planar_constraint_along_direction, plane_incidence, rail_incidence, roller_flat, Line, Motor,
    Plane, Point, SupportFlatKind,
};

const TOL: f64 = 1e-10;

fn approx(a: f64, b: f64) -> bool {
    (a - b).abs() <= TOL
}

fn xyz_close(a: [f64; 3], b: [f64; 3]) -> bool {
    a.iter().zip(b).all(|(u, v)| approx(*u, v))
}

#[test]
fn join_of_two_points_is_the_line_through_them() {
    let a = Point::euclidean(0.0, 0.0, 1.0);
    let b = Point::euclidean(2.0, 0.0, 1.0);
    let line = a.join(b);
    assert!(approx(line.distance(Point::euclidean(1.0, 0.0, 1.0)), 0.0));
    assert!(line.distance(Point::euclidean(1.0, 1.0, 1.0)) > 0.5);
    let d = line.direction;
    let m = line.moment;
    assert!(approx(d[0] * m[0] + d[1] * m[1] + d[2] * m[2], 0.0));
}

#[test]
fn meet_of_two_planes_is_their_intersection_line() {
    let xy = Plane::new(0.0, 0.0, 1.0, 0.0); // z = 0
    let xz = Plane::new(0.0, 1.0, 0.0, 0.0); // y = 0
    let line = xy.meet(xz);
    assert!(approx(line.distance(Point::euclidean(3.0, 0.0, 0.0)), 0.0));
    assert!(line.distance(Point::euclidean(0.0, 1.0, 0.0)) > 0.5);
}

#[test]
fn three_planes_meet_at_a_point() {
    let x0 = Plane::new(1.0, 0.0, 0.0, -2.0); // x = 2
    let y0 = Plane::new(0.0, 1.0, 0.0, -3.0); // y = 3
    let z0 = Plane::new(0.0, 0.0, 1.0, -4.0); // z = 4
    let p = x0.meet_line(y0.meet(z0)).xyz().unwrap();
    assert!(xyz_close(p, [2.0, 3.0, 4.0]));
}

#[test]
fn three_points_span_a_plane() {
    let plane = Plane::through(
        Point::euclidean(0.0, 0.0, 0.0),
        Point::euclidean(1.0, 0.0, 0.0),
        Point::euclidean(0.0, 1.0, 0.0),
    );
    assert!(approx(
        plane.signed_distance(Point::euclidean(0.3, 0.4, 0.0)),
        0.0
    ));
    assert!(approx(
        plane.signed_distance(Point::euclidean(0.0, 0.0, 2.0)).abs(),
        2.0
    ));
}

#[test]
fn translator_moves_the_origin() {
    let t = Motor::translator(1.5, -2.0, 4.0);
    let p = t
        .apply_point(Point::euclidean(0.0, 0.0, 0.0))
        .xyz()
        .unwrap();
    assert!(xyz_close(p, [1.5, -2.0, 4.0]));
}

#[test]
fn rotator_around_z_sends_x_to_y() {
    let r = Motor::rotator([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);
    let p = r
        .apply_point(Point::euclidean(1.0, 0.0, 0.0))
        .xyz()
        .unwrap();
    assert!(xyz_close(p, [0.0, 1.0, 0.0]));
}

#[test]
fn rail_incidence_vanishes_on_the_supporting_line() {
    let start = [0.0, 0.0, 0.0];
    let end = [4.0, 2.0, 0.0];
    assert!(approx(rail_incidence([2.0, 1.0, 0.0], start, end), 0.0));
    assert!(rail_incidence([2.0, 1.0, 1.0], start, end) > 0.5);
}

#[test]
fn plane_incidence_matches_signed_distance() {
    let origin = [0.0, 0.0, 3.0];
    let x = [1.0, 0.0, 0.0];
    let y = [0.0, 1.0, 0.0];
    assert!(approx(plane_incidence([1.0, 2.0, 3.0], origin, x, y), 0.0));
    assert!(approx(plane_incidence([0.0, 0.0, 5.0], origin, x, y), 2.0));
}

#[test]
fn planar_constraint_along_normal_is_scaled_incidence() {
    let origin = [0.0, 0.0, 0.0];
    let x = [1.0, 0.0, 0.0];
    let y = [0.0, 1.0, 0.0];
    let n = [0.0, 0.0, 2.0];
    let node = [1.0, 1.0, 4.0];
    let t = planar_constraint_along_direction(node, origin, x, y, n).unwrap();
    let incidence = plane_incidence(node, origin, x, y);
    // t = n·(O-P)/(n·n_dir). With plane z=0 and d=(0,0,2): t = (0-4)/2 = -2
    // Euclidean incidence is +4. They agree only up to a chart, not as values.
    assert!(approx(t, -2.0));
    assert!(approx(incidence, 4.0));
}

#[test]
fn roller_axis_counts_map_to_flats() {
    assert_eq!(roller_flat([false, false, false]), SupportFlatKind::Point);
    assert_eq!(roller_flat([true, false, false]), SupportFlatKind::Line);
    assert_eq!(roller_flat([true, true, false]), SupportFlatKind::Plane);
    assert_eq!(roller_flat([true, true, true]), SupportFlatKind::NotAFlat);
}

#[test]
fn line_parameter_is_the_rail_chart() {
    let start = [1.0, 1.0, 1.0];
    let end = [5.0, 1.0, 1.0];
    let line = Line::from_segment(start, end);
    let t = line
        .segment_parameter(Point::euclidean(3.0, 1.0, 1.0), start, end)
        .unwrap();
    assert!(approx(t, 0.5));
}
