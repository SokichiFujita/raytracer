use std::f32::INFINITY;

use crate::material::Material;
use na::Matrix4;
use ulid::Ulid;

#[derive(Clone, Debug)]
pub struct Cylinder {
    pub id: String,
    pub transformation: Matrix4<f32>,
    pub material: Material,
    pub min: f32,
    pub max: f32,
    pub closed: bool,
}

impl PartialEq for Cylinder {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Cylinder {
    pub fn new(
        transformation: Option<Matrix4<f32>>,
        matrial: Option<Material>,
        min: Option<f32>,
        max: Option<f32>,
        closed: Option<bool>,
    ) -> Cylinder {
        Cylinder {
            id: Ulid::new().to_string(),
            transformation: match transformation {
                Some(x) => x,
                None => Matrix4::<f32>::identity(),
            },
            material: match matrial {
                Some(x) => x,
                None => Material::new_default(),
            },
            min: match min {
                Some(x) => x,
                None => -INFINITY,
            },
            max: match max {
                Some(x) => x,
                None => INFINITY,
            },
            closed: match closed {
                Some(x) => x,
                None => false,
            },
        }
    }
    pub fn new_default() -> Cylinder {
        Cylinder::new(None, None, None, None, None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{epsilon::EPSILON, ray::Ray, shape::Shape, tuple::TupleOperation};
    use approx::assert_relative_eq;
    use nalgebra::Vector4;

    #[test]
    fn ray_misses_cylinder() {
        let cylinder = Shape::Cylinder(Cylinder::new_default());
        let test_set = vec![
            (Vector4::point(1., 0., 0.), Vector4::vector(0., 1., 0.)),
            (Vector4::point(0., 0., 0.), Vector4::vector(0., 1., 0.)),
            (Vector4::point(0., 0., -5.), Vector4::vector(1., 1., 1.)),
        ];
        for (origin, direction) in test_set {
            let ray = Ray::new(origin, direction.normalize());
            let xs = cylinder.intersect(&ray);
            assert_eq!(xs.len(), 0);
        }
    }

    #[test]
    fn ray_intersects_cylinder() {
        let cylinder = Shape::Cylinder(Cylinder::new_default());
        let test_set = vec![
            (
                Vector4::point(1., 0., -5.),
                Vector4::vector(0., 0., 1.),
                5.,
                5.,
            ),
            (
                Vector4::point(0., 0., -5.),
                Vector4::vector(0., 0., 1.),
                4.,
                6.,
            ),
            (
                Vector4::point(0.5, 0., -5.),
                Vector4::vector(0.1, 1., 1.),
                6.80798,
                7.08872,
            ),
        ];
        for (origin, direction, t0, t1) in test_set {
            let ray = Ray::new(origin, direction.normalize());
            let xs = cylinder.intersect(&ray);
            assert_eq!(xs.len(), 2);
            assert_relative_eq!(xs[0].t, t0, epsilon = EPSILON);
            assert_relative_eq!(xs[1].t, t1, epsilon = EPSILON);
        }
    }

    #[test]
    fn intersecting_constrained_cylinder() {
        let cylinder = Shape::Cylinder(Cylinder::new_default());
        let test_set = vec![
            (
                Vector4::point(1., 0., -5.),
                Vector4::vector(0., 0., 1.),
                5.,
                5.,
            ),
            (
                Vector4::point(0., 0., -5.),
                Vector4::vector(0., 0., 1.),
                4.,
                6.,
            ),
            (
                Vector4::point(0.5, 0., -5.),
                Vector4::vector(0.1, 1., 1.),
                6.80798,
                7.08872,
            ),
        ];
        for (origin, direction, t0, t1) in test_set {
            let ray = Ray::new(origin, direction.normalize());
            let xs = cylinder.intersect(&ray);
            assert_eq!(xs.len(), 2);
            assert_relative_eq!(xs[0].t, t0, epsilon = EPSILON);
            assert_relative_eq!(xs[1].t, t1, epsilon = EPSILON);
        }
    }

    #[test]
    fn normal_cylinder() {
        for (point, normal) in vec![
            (Vector4::point(1., 0., 0.), Vector4::vector(1., 0., 0.)),
            (Vector4::point(0., 5., -1.), Vector4::vector(0., 0., -1.)),
            (Vector4::point(0., 2., -1.), Vector4::vector(0., 0., -1.)),
            (Vector4::point(-1., 1., 0.), Vector4::vector(-1., 0., 0.)),
        ] {
            let shape = Shape::Cylinder(Cylinder::new_default());
            let n = shape.normal(point);
            assert_relative_eq!(n, normal);
        }
    }

    #[test]
    fn intersect_constrained_cylinder() {
        for (origin, direction, count) in vec![
            (Vector4::point(0., 1.5, 0.), Vector4::vector(0.1, 1., 0.), 0),
            (Vector4::point(0., 3., -5.), Vector4::vector(0., 0., 1.), 0),
            (Vector4::point(0., 0., -5.), Vector4::vector(0., 0., 1.), 0),
            (Vector4::point(0., 2., -5.), Vector4::vector(0., 0., 1.), 0),
            (Vector4::point(0., 1., -5.), Vector4::vector(0., 0., 1.), 0),
            (Vector4::point(0., 1.5, -2.), Vector4::vector(0., 0., 1.), 2),
        ] {
            let shape =
                Shape::Cylinder(Cylinder::new(None, None, Some(1.0), Some(2.0), Some(false)));
            let ray = Ray::new(origin, direction.normalize());
            let xs = shape.intersect(&ray);
            assert_eq!(xs.len(), count);
        }
    }

    #[test]
    fn intersect_caps_closed_cylinder() {
        for (origin, direction, count) in vec![
            (Vector4::point(0., 3., 0.), Vector4::vector(0., -1., 0.), 2),
            (Vector4::point(0., 3., -2.), Vector4::vector(0., -1., 2.), 2),
            (Vector4::point(0., 4., -2.), Vector4::vector(0., -1., 1.), 2),
            (Vector4::point(0., 0., -2.), Vector4::vector(0., 1., 2.), 2),
            (Vector4::point(0., -1., -2.), Vector4::vector(0., 1., 1.), 2),
        ] {
            let shape =
                Shape::Cylinder(Cylinder::new(None, None, Some(1.0), Some(2.0), Some(true)));
            let ray = Ray::new(origin, direction.normalize());
            let xs = shape.intersect(&ray);
            assert_eq!(xs.len(), count);
        }
    }

    #[test]
    fn normal_closed() {
        for (point, normal) in vec![
            (Vector4::point(0., 1., 0.), Vector4::vector(0., -1., 0.)),
            (Vector4::point(0.5, 1., 0.), Vector4::vector(0., -1., 0.)),
            (Vector4::point(0., 1., 0.5), Vector4::vector(0., -1., 0.)),
            (Vector4::point(0., 2., 0.), Vector4::vector(0., 1., 0.)),
            (Vector4::point(0.5, 2., 0.), Vector4::vector(0., 1., 0.)),
            (Vector4::point(0., 2., 0.5), Vector4::vector(0., 1., 0.)),
        ] {
            let shape =
                Shape::Cylinder(Cylinder::new(None, None, Some(1.0), Some(2.0), Some(true)));
            assert_relative_eq!(shape.normal(point), normal);
        }
    }
}
