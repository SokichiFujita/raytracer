use na::{Matrix4, Vector4};

use crate::{material::Material, matrix::CGMatrix, shape::Shape, tuple::TupleOperation};

#[derive(Clone, Debug)]
pub struct Triangle {
    pub id: String,
    pub transformation: Matrix4<f32>,
    pub material: Material,
    pub p1: Vector4<f32>, //top
    pub p2: Vector4<f32>, //left bottom
    pub p3: Vector4<f32>, //right bottom
    pub e1: Vector4<f32>, //p2 - p1 left
    pub e2: Vector4<f32>, //p3 - p1 right
    pub normal: Vector4<f32>,
}

impl PartialEq for Triangle {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Triangle {
    pub fn new(
        transformation: Option<Matrix4<f32>>,
        matrial: Option<Material>,
        p1: Vector4<f32>,
        p2: Vector4<f32>,
        p3: Vector4<f32>,
    ) -> Triangle {
        let e1 = p2 - p1;
        let e2 = p3 - p2;
        Triangle {
            id: Shape::generate_id(Some("cube")),
            transformation: match transformation {
                Some(x) => x,
                None => Matrix4::<f32>::identity(),
            },
            material: match matrial {
                Some(x) => x,
                None => Material::new_default(),
            },
            p1,
            p2,
            p3,
            e1: p2 - p1,
            e2: p3 - p1,
            normal: e2.cross4(&e1).normalize(),
        }
    }
    pub fn new_default() -> Triangle {
        Triangle::new(
            None,
            None,
            Vector4::point(0., 1., 0.),
            Vector4::point(-1., 0., 0.),
            Vector4::point(1., 0., 0.),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ray::Ray, shape::Shape, tuple::TupleOperation};
    use approx::assert_relative_eq;
    use nalgebra::Vector4;
    #[test]
    fn constructing_triangle() {
        let p1 = Vector4::point(0., 1., 0.);
        let p2 = Vector4::point(-1., 0., 0.);
        let p3 = Vector4::point(1., 0., 0.);
        let t = Triangle::new(None, None, p1, p2, p3);
        assert_eq!(t.p1, p1);
        assert_eq!(t.p2, p2);
        assert_eq!(t.p3, p3);
        assert_eq!(t.e1, Vector4::vector(-1., -1., 0.));
        assert_eq!(t.e2, Vector4::vector(1., -1., 0.));
        assert_eq!(t.normal, Vector4::vector(0., 0., -1.));
    }

    #[test]
    fn finding_normal_on_triangle() {
        let p1 = Vector4::point(0., 1., 0.);
        let p2 = Vector4::point(-1., 0., 0.);
        let p3 = Vector4::point(1., 0., 0.);
        let t = Shape::Triangle(Triangle::new(None, None, p1, p2, p3));
        let n1 = t.local_normal(Vector4::point(0., 0.5, 0.));
        let n2 = t.local_normal(Vector4::point(-0.5, 0.75, 0.));
        let n3 = t.local_normal(Vector4::point(0.5, 0.25, 0.));

        assert_eq!(t.triangle().unwrap().normal, n1);
        assert_eq!(t.triangle().unwrap().normal, n2);
        assert_eq!(t.triangle().unwrap().normal, n3);
    }

    #[test]
    fn ray_misses_p1_pe_edge() {
        let t = Shape::Triangle(Triangle::new(
            None,
            None,
            Vector4::point(0., 1., 0.),
            Vector4::point(-1., 0., 0.),
            Vector4::point(1., 0., 0.),
        ));
        let ray = Ray::new(Vector4::point(1., 1., -2.), Vector4::point(0., 0., 1.));
        let xs = t.local_intersect(&ray);
        assert_eq!(xs.len(), 0);
    }

    #[test]
    fn ray_misses_p1_p3_edge() {
        let t = Shape::Triangle(Triangle::new(
            None,
            None,
            Vector4::point(0., 1., 0.),
            Vector4::point(-1., 0., 0.),
            Vector4::point(1., 0., 0.),
        ));
        let ray = Ray::new(Vector4::point(1., 1., -2.), Vector4::point(0., 0., 1.));
        let xs = t.local_intersect(&ray);
        assert_eq!(xs.len(), 0);
    }

    #[test]
    fn ray_misses_p1_p2_edge() {
        let t = Shape::Triangle(Triangle::new(
            None,
            None,
            Vector4::point(0., 1., 0.),
            Vector4::point(-1., 0., 0.),
            Vector4::point(1., 0., 0.),
        ));
        let ray = Ray::new(Vector4::point(-1., 1., -2.), Vector4::point(0., 0., 1.));
        let xs = t.local_intersect(&ray);
        assert_eq!(xs.len(), 0);
    }

    #[test]
    fn ray_misses_p2_p3_edge() {
        let t = Shape::Triangle(Triangle::new(
            None,
            None,
            Vector4::point(0., 1., 0.),
            Vector4::point(-1., 0., 0.),
            Vector4::point(1., 0., 0.),
        ));
        let ray = Ray::new(Vector4::point(0., -1., -2.), Vector4::point(0., 0., 1.));
        let xs = t.local_intersect(&ray);
        assert_eq!(xs.len(), 0);
    }

    #[test]
    fn ray_strikes_triangle() {
        let t = Shape::Triangle(Triangle::new(
            None,
            None,
            Vector4::point(0., 1., 0.),
            Vector4::point(-1., 0., 0.),
            Vector4::point(1., 0., 0.),
        ));
        let ray = Ray::new(Vector4::point(0., 0.5, -2.), Vector4::point(0., 0., 1.));
        let xs = t.local_intersect(&ray);
        assert_eq!(xs.len(), 1);
        assert_eq!(xs[0].t, 2.0);
    }
}
