use na::Matrix4;

use crate::{material::Material, shape::Shape};

#[derive(Clone, Debug)]
pub struct Cube {
    pub id: String,
    pub transformation: Matrix4<f32>,
    pub material: Material,
}

impl PartialEq for Cube {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Cube {
    pub fn new(transformation: Option<Matrix4<f32>>, matrial: Option<Material>) -> Cube {
        Cube {
            id: Shape::generate_id(Some("cube")),
            transformation: match transformation {
                Some(x) => x,
                None => Matrix4::<f32>::identity(),
            },
            material: match matrial {
                Some(x) => x,
                None => Material::new_default(),
            },
        }
    }
    pub fn new_default() -> Cube {
        Cube::new(None, None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ray::Ray, shape::Shape, tuple::TupleOperation};
    use approx::assert_relative_eq;
    use nalgebra::Vector4;
    #[test]
    fn ray_intersects_cube() {
        let cube = Shape::Cube(Cube::new_default());
        let test_set = vec![
            (
                Vector4::point(5., 0.5, 0.),
                Vector4::vector(-1., 0., 0.),
                4.,
                6.,
            ),
            (
                Vector4::point(-5., 0.5, 0.),
                Vector4::vector(1., 0., 0.),
                4.,
                6.,
            ),
            (
                Vector4::point(0.5, 5., 0.),
                Vector4::vector(0., -1., 0.),
                4.,
                6.,
            ),
            (
                Vector4::point(0.5, -5., 0.),
                Vector4::vector(0., 1., 0.),
                4.,
                6.,
            ),
            (
                Vector4::point(0.5, 0., 5.),
                Vector4::vector(0., 0., -1.),
                4.,
                6.,
            ),
            (
                Vector4::point(0.5, 0., -5.),
                Vector4::vector(0., 0., 1.),
                4.,
                6.,
            ),
            (
                Vector4::point(0., 0.5, 0.),
                Vector4::vector(0., 0., 1.),
                -1.,
                1.,
            ),
        ];
        for (origin, direction, t0, t1) in test_set {
            let ray = Ray::new(origin, direction);
            let xs = cube.intersect(&ray);
            assert_relative_eq!(xs[0].t, t0);
            assert_relative_eq!(xs[1].t, t1);
        }
    }

    #[test]
    fn ray_misses_cube() {
        let cube = Shape::Cube(Cube::new_default());
        let test_set = vec![
            (
                Vector4::point(-2., 0., 0.),
                Vector4::vector(0.2673, 0.5345, 0.8018),
            ),
            (
                Vector4::point(0., -2., 0.),
                Vector4::vector(0.8018, 0.2673, 0.5345),
            ),
            (
                Vector4::point(0., 0., -2.),
                Vector4::vector(0.5345, 0.8018, 0.2673),
            ),
            (Vector4::point(2., 0., 2.), Vector4::vector(0., 0., -1.)),
            (Vector4::point(0., 2., 2.), Vector4::vector(0., -1., 0.)),
            (Vector4::point(2., 2., 0.), Vector4::vector(-1., 0., 0.)),
        ];
        for (origin, direction) in test_set {
            let ray = Ray::new(origin, direction);
            let xs = cube.intersect(&ray);
            assert_eq!(xs.len(), 0);
        }
    }

    #[test]
    fn normal_on_surface_of_cube() {
        let cube = Shape::Cube(Cube::new_default());
        let test_set = vec![
            (Vector4::point(1., 0.5, -0.8), Vector4::vector(1., 0., 0.)),
            (Vector4::point(-1., -0.2, 0.9), Vector4::vector(-1., 0., 0.)),
            (Vector4::point(-0.4, 1., -0.1), Vector4::vector(0., 1., 0.)),
            (Vector4::point(0.3, -1., -0.7), Vector4::vector(0., -1., 0.)),
            (Vector4::point(-0.6, 0.3, 1.), Vector4::vector(0., 0., 1.)),
            (Vector4::point(0.4, 0.4, -1.), Vector4::vector(0., 0., -1.)),
            (Vector4::point(1., 1., 1.), Vector4::vector(1., 0., 0.)),
            (Vector4::point(-1., -1., -1.), Vector4::vector(-1., 0., 0.)),
        ];
        for (point, normal) in test_set {
            assert_eq!(cube.normal(point), normal);
        }
    }
}
