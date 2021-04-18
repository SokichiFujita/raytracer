use na::{Matrix4, Vector4};

use crate::{
    epsilon::EPSILON, intersection::Intersection, material::Material, matrix::CGMatrix,
    shape::Shape, tuple::TupleOperation,
};

#[derive(Clone, Debug)]
pub struct SmoothTriangle {
    pub id: String,
    pub transformation: Matrix4<f32>,
    pub material: Material,
    pub p1: Vector4<f32>, //top
    pub p2: Vector4<f32>, //left bottom
    pub p3: Vector4<f32>, //right bottom
    pub n1: Vector4<f32>, //normal at top
    pub n2: Vector4<f32>, //normal at left bottom
    pub n3: Vector4<f32>, //normal at right bottom
    pub e1: Vector4<f32>, //p2 - p1 left
    pub e2: Vector4<f32>, //p3 - p1 right
}

impl PartialEq for SmoothTriangle {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl SmoothTriangle {
    pub fn new(
        transformation: Option<Matrix4<f32>>,
        matrial: Option<Material>,
        p1: Vector4<f32>,
        p2: Vector4<f32>,
        p3: Vector4<f32>,
    ) -> SmoothTriangle {
        SmoothTriangle::new_with_n(
            transformation,
            matrial,
            p1,
            p2,
            p3,
            Vector4::vector(0., 0., 0.),
            Vector4::vector(0., 0., 0.),
            Vector4::vector(0., 0., 0.),
        )
    }
    pub fn new_default() -> SmoothTriangle {
        SmoothTriangle::new(
            None,
            None,
            Vector4::point(0., 1., 0.),
            Vector4::point(-1., 0., 0.),
            Vector4::point(1., 0., 0.),
        )
    }
    pub fn new_with_n(
        transformation: Option<Matrix4<f32>>,
        matrial: Option<Material>,
        p1: Vector4<f32>,
        p2: Vector4<f32>,
        p3: Vector4<f32>,
        n1: Vector4<f32>,
        n2: Vector4<f32>,
        n3: Vector4<f32>,
    ) -> SmoothTriangle {
        let e1 = p2 - p1;
        let e2 = p3 - p2;
        SmoothTriangle {
            id: Shape::generate_id(Some("smoothtriangle")),
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
            n1,
            n2,
            n3,
            e1: p2 - p1,
            e2: p3 - p1,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ray::Ray, shape::Shape, tuple::TupleOperation};
    use approx::assert_relative_eq;
    use nalgebra::Vector4;
    #[test]
    fn constructing_smooth_triangle() {
        let p1 = Vector4::point(0., 1., 0.);
        let p2 = Vector4::point(-1., 0., 0.);
        let p3 = Vector4::point(1., 0., 0.);
        let n1 = Vector4::vector(0., 1., 0.);
        let n2 = Vector4::vector(-1., 0., 0.);
        let n3 = Vector4::vector(1., 0., 0.);

        let tri = Shape::SmoothTriangle(SmoothTriangle::new_with_n(
            None, None, p1, p2, p3, n1, n2, n3,
        ));
        assert_eq!(tri.smooth_triangle().unwrap().p1, p1);
        assert_eq!(tri.smooth_triangle().unwrap().p2, p2);
        assert_eq!(tri.smooth_triangle().unwrap().p3, p3);
        assert_eq!(tri.smooth_triangle().unwrap().n1, n1);
        assert_eq!(tri.smooth_triangle().unwrap().n2, n2);
        assert_eq!(tri.smooth_triangle().unwrap().n3, n3);
    }

    #[test]
    fn intersection_can_encapsulate_u_v() {
        let p1 = Vector4::point(0., 1., 0.);
        let p2 = Vector4::point(-1., 0., 0.);
        let p3 = Vector4::point(1., 0., 0.);
        let n1 = Vector4::vector(0., 1., 0.);
        let n2 = Vector4::vector(-1., 0., 0.);
        let n3 = Vector4::vector(1., 0., 0.);

        let tri = Shape::SmoothTriangle(SmoothTriangle::new_with_n(
            None, None, p1, p2, p3, n1, n2, n3,
        ));
        let i = Intersection::new_with_uv(3.5, &tri, 0.2, 0.4);
        assert_eq!(i.u.unwrap(), 0.2);
        assert_eq!(i.v.unwrap(), 0.4);
    }

    #[test]
    fn intersection_with_smooth_triangle_stores_u_v() {
        let p1 = Vector4::point(0., 1., 0.);
        let p2 = Vector4::point(-1., 0., 0.);
        let p3 = Vector4::point(1., 0., 0.);
        let n1 = Vector4::vector(0., 1., 0.);
        let n2 = Vector4::vector(-1., 0., 0.);
        let n3 = Vector4::vector(1., 0., 0.);

        let tri = Shape::SmoothTriangle(SmoothTriangle::new_with_n(
            None, None, p1, p2, p3, n1, n2, n3,
        ));
        let ray = Ray::new(Vector4::point(-0.2, 0.3, -2.0), Vector4::vector(0., 0., 1.));
        let xs = tri.local_intersect(&ray);
        assert_eq!(xs[0].u.unwrap(), 0.45);
        assert_eq!(xs[0].v.unwrap(), 0.25);
    }

    #[test]
    fn smooth_triangle_uses_u_v_to_interpolate_normal() {
        let p1 = Vector4::point(0., 1., 0.);
        let p2 = Vector4::point(-1., 0., 0.);
        let p3 = Vector4::point(1., 0., 0.);
        let n1 = Vector4::vector(0., 1., 0.);
        let n2 = Vector4::vector(-1., 0., 0.);
        let n3 = Vector4::vector(1., 0., 0.);
        let tri = Shape::SmoothTriangle(SmoothTriangle::new_with_n(
            None, None, p1, p2, p3, n1, n2, n3,
        ));
        let i = Intersection::new_with_uv(1., &tri, 0.45, 0.25);
        let n = tri.normal(Vector4::point(0., 0., 0.), Some(i));
        assert_relative_eq!(n, Vector4::vector(-0.5547, 0.83205, 0.0), epsilon = EPSILON);
    }
}
