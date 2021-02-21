use std::f32::INFINITY;

use crate::{material::Material, shape::Shape};
use na::Matrix4;

#[derive(Clone, Debug)]
pub struct Cone {
    pub id: String,
    pub transformation: Matrix4<f32>,
    pub material: Material,
    pub min: f32,
    pub max: f32,
    pub closed: bool,
    pub parent: Option<String>,
}

impl PartialEq for Cone {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Cone {
    pub fn new(
        transformation: Option<Matrix4<f32>>,
        matrial: Option<Material>,
        min: Option<f32>,
        max: Option<f32>,
        closed: Option<bool>,
        parent: Option<String>,
    ) -> Cone {
        Cone {
            id: Shape::generate_id(Some("cone")),
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
            parent,
        }
    }
    pub fn new_default() -> Cone {
        Cone::new(None, None, None, None, None, None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{epsilon::EPSILON, ray::Ray, shape::Shape, tuple::TupleOperation};
    use approx::assert_relative_eq;
    use nalgebra::Vector4;

    #[test]
    fn intersecting_cone_ray() {
        let test_sets = vec![
            (
                Vector4::point(0.0, 0.0, -5.0),
                Vector4::vector(0.0, 0.0, 1.0),
                5.0,
                5.0,
            ),
            (
                Vector4::point(0.0, 0.0, -5.0),
                Vector4::vector(1.0, 1.0, 1.0),
                8.66025,
                8.66025,
            ),
            (
                Vector4::point(1.0, 1.0, -5.0),
                Vector4::vector(-0.5, -1.0, 1.0),
                4.55006,
                49.44994,
            ),
        ];
        let shape = Shape::Cone(Cone::new_default());
        for (origin, direction, t0, t1) in test_sets.into_iter() {
            let xs = shape.intersect(&Ray::new(origin, direction.normalize()));
            assert_eq!(xs.len(), 2);
            assert_relative_eq!(xs[0].t, t0, epsilon = EPSILON);
            assert_relative_eq!(xs[1].t, t1, epsilon = EPSILON);
        }
    }

    #[test]
    fn intersecting_cone_ray_parallel_to_halves() {
        let test_sets = vec![(
            Vector4::point(0., 0., -1.),
            Vector4::vector(0., 1., 1.),
            0.35355,
        )];
        let shape = Shape::Cone(Cone::new_default());
        for (origin, direction, t0) in test_sets.into_iter() {
            let xs = shape.intersect(&Ray::new(origin, direction.normalize()));
            assert_eq!(xs.len(), 1);
            assert_relative_eq!(xs[0].t, t0, epsilon = EPSILON);
        }
    }

    #[test]
    fn intersecting_cone_end_caps() {
        let test_sets = vec![
            (
                Vector4::point(0.0, 0.0, -5.0),
                Vector4::vector(0.0, 1.0, 0.0),
                0,
            ),
            (
                Vector4::point(0.0, 0.0, -0.25),
                Vector4::vector(0.0, 1.0, 1.0),
                2,
            ),
            (
                Vector4::point(0.0, 0.0, -0.25),
                Vector4::vector(0.0, 1.0, 0.0),
                4,
            ),
        ];
        let shape = Shape::Cone(Cone::new(
            None,
            None,
            Some(-0.5),
            Some(0.5),
            Some(true),
            None,
        ));
        for (origin, direction, count) in test_sets.into_iter() {
            let xs = shape.intersect(&Ray::new(origin, direction.normalize()));
            assert_eq!(xs.len(), count);
        }
    }

    #[test]
    fn normal_vector_cone() {
        let test_sets = vec![
            (
                Vector4::point(0.0, 0.0, 0.0),
                Vector4::vector(0.0, 0.0, 0.0),
            ),
            (
                Vector4::point(1.0, 1.0, 1.0),
                Vector4::vector(1.0, -(2.0_f32.sqrt()), 1.0),
            ),
            (
                Vector4::point(-1.0, -1.0, 0.0),
                Vector4::vector(-1.0, 1.0, 0.0),
            ),
        ];
        let shape = Shape::Cone(Cone::new_default());
        for (point, normal) in test_sets.into_iter() {
            assert_relative_eq!(shape.normal(point), normal);
        }
    }
}
