use crate::{material::Material, shape::Shape};
use na::{Matrix4, Vector4};

use crate::tuple::TupleOperation;

#[derive(Clone, Debug)]
pub struct Sphere {
    pub id: String,
    pub transformation: Matrix4<f32>,
    pub material: Material,
    pub position: Vector4<f32>, //origin
    pub r: f32,
    pub parent: Option<String>,
}

impl PartialEq for Sphere {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Sphere {
    pub fn new(
        transformation: Option<Matrix4<f32>>,
        matrial: Option<Material>,
        position: Option<Vector4<f32>>,
        r: Option<f32>,
        parent: Option<String>,
    ) -> Sphere {
        Sphere {
            id: Shape::generate_id(Some("sphere")),
            transformation: match transformation {
                Some(x) => x,
                None => Matrix4::<f32>::identity(),
            },
            material: match matrial {
                Some(x) => x,
                None => Material::new_default(),
            }, //origin
            position: match position {
                Some(x) => x,
                None => Vector4::point(0.0, 0.0, 0.0),
            }, //origin
            r: match r {
                Some(x) => x,
                None => 1.0,
            },
            parent,
        }
    }
    pub fn new_default() -> Sphere {
        Sphere::new(None, None, None, None, None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{epsilon::EPSILON, matrix::CGMatrix, ray::Ray, shape::Shape};
    use approx::assert_relative_eq;
    use std::f32::consts::PI;

    //Ch5 P.59
    #[test]
    fn intersect_sphere_two_points() {
        let sphere = Shape::Sphere(Sphere::new_default());
        let ray = Ray::new(
            Vector4::point(0.0, 0.0, -5.0),
            Vector4::vector(0.0, 0.0, 1.0),
        );
        let intersections = sphere.intersect(&ray);
        let length = intersections.iter().count();
        assert_eq!(length, 2)
    }

    //Ch5 P.60
    #[test]
    fn intersect_sphere_tangent() {
        let sphere = Shape::Sphere(Sphere::new_default());
        let ray = Ray::new(
            Vector4::point(0.0, 1.0, -5.0),
            Vector4::vector(0.0, 0.0, 1.0),
        );
        let intersections = sphere.intersect(&ray);
        let length = intersections.iter().count();
        assert_eq!(length, 2)
    }

    //Ch5 P.59
    #[test]
    fn intersect_sphere_no_point() {
        let sphere = Shape::Sphere(Sphere::new_default());
        let ray = Ray::new(
            Vector4::point(0.0, 2.0, -5.0),
            Vector4::vector(0.0, 0.0, 1.0),
        );
        let intersections = sphere.intersect(&ray);
        let length = intersections.iter().count();
        assert_eq!(length, 0)
    }

    //Ch5 P.61
    #[test]
    fn originates_inside_sphere() {
        let sphere = Shape::Sphere(Sphere::new_default());
        let ray = Ray::new(
            Vector4::point(0.0, 0.0, 0.0),
            Vector4::vector(0.0, 0.0, 1.0),
        );
        let intersections = sphere.intersect(&ray);
        let length = intersections.iter().count();
        assert_eq!(length, 2)
    }

    //Ch5 P.69
    #[test]
    fn sphere_default() {
        let sphere = Sphere::new_default();
        let identity = Matrix4::<f32>::identity();
        assert_eq!(sphere.transformation, identity)
    }

    //Ch5 P.69
    #[test]
    fn change_sphere_transformation() {
        let mut sphere = Sphere::new_default();
        let translation = Matrix4::<f32>::translation(2.0, 3.0, 4.0);
        sphere.transformation = translation;
        assert_eq!(sphere.transformation, translation)
    }

    //Ch5 P.69
    #[test]
    fn intersecting_scaled_sphere_with_ray() {
        let ray = Ray::from_tuple((0.0, 0.0, -5.0), (0.0, 0.0, 1.0));
        let mut sphere = Sphere::new_default();
        let scaling = Matrix4::scaling(2.0, 2.0, 2.0);
        sphere.transformation = scaling;
        assert_eq!(sphere.transformation, scaling);

        let shape = Shape::Sphere(sphere);
        let xs = shape.intersect(&ray);
        assert_eq!(xs.len(), 2);
        assert_eq!(xs.get(0).unwrap().t, 3.0);
        assert_eq!(xs.get(1).unwrap().t, 7.0);
    }

    //Ch5 P.69
    #[test]
    fn intersecting_translated_sphere_with_ray() {
        let ray = Ray::from_tuple((0.0, 0.0, -5.0), (0.0, 0.0, 1.0));
        let mut sphere = Sphere::new_default();
        let translation = Matrix4::translation(5.0, 0.0, 0.0);
        sphere.transformation = translation;
        assert_eq!(sphere.transformation, translation);

        let shape = Shape::Sphere(sphere);
        let xs = shape.intersect(&ray);
        assert_eq!(xs.len(), 0);
    }

    //Ch6 P.78
    #[test]
    fn normal_on_sphere_at_xaxis() {
        let sphere = Shape::Sphere(Sphere::new_default());
        let point = Vector4::point(1., 0., 0.);
        let normal = sphere.normal(point);
        assert_eq!(normal, Vector4::vector(1., 0., 0.));
    }

    //Ch6 P.78
    #[test]
    fn normal_on_sphere_at_yaxis() {
        let sphere = Shape::Sphere(Sphere::new_default());
        let point = Vector4::point(0., 1., 0.);
        let normal = sphere.normal(point);
        assert_eq!(normal, Vector4::vector(0., 1., 0.));
    }

    //Ch6 P.78
    #[test]
    fn normal_on_sphere_at_zaxis() {
        let sphere = Shape::Sphere(Sphere::new_default());
        let point = Vector4::point(0., 0., 1.);
        let normal = sphere.normal(point);
        assert_eq!(normal, Vector4::vector(0., 0., 1.));
    }

    //Ch6 P.78
    #[test]
    fn normal_on_sphere_at_non_axis() {
        let sphere = Shape::Sphere(Sphere::new_default());
        let point = Vector4::point(
            3.0_f32.sqrt() / 3.0,
            3.0_f32.sqrt() / 3.0,
            3.0_f32.sqrt() / 3.0,
        );
        let normal = sphere.normal(point);
        assert_relative_eq!(
            normal,
            Vector4::vector(
                3.0_f32.sqrt() / 3.0,
                3.0_f32.sqrt() / 3.0,
                3.0_f32.sqrt() / 3.0
            )
        );
    }

    //Ch6 P.78
    #[test]
    fn normal_is_normalized() {
        let sphere = Shape::Sphere(Sphere::new_default());
        let point = Vector4::point(
            3.0_f32.sqrt() / 3.0,
            3.0_f32.sqrt() / 3.0,
            3.0_f32.sqrt() / 3.0,
        );
        let normal = sphere.normal(point);
        assert_relative_eq!(normal, normal.normalize());
    }

    //Ch6 P.80
    #[test]
    fn computing_normal_on_translated_sphere() {
        let sphere = Shape::Sphere(Sphere::new(
            Some(Matrix4::translation(0., 1., 0.)),
            None,
            None,
            None,
            None,
        ));

        let point = Vector4::point(0., 1.70711, -0.70711);
        let normal = sphere.normal(point);
        assert_relative_eq!(
            normal,
            Vector4::vector(0., 0.70711, -0.70711),
            epsilon = EPSILON
        );
    }

    //Ch6 P.80
    #[test]
    fn computing_normal_on_transformed_sphere() {
        let sphere = Shape::Sphere(Sphere::new(
            Some(Matrix4::scaling(1.0, 0.5, 1.0) * Matrix4::rotation_z(PI / 5.0)),
            None,
            None,
            None,
            None,
        ));

        let point = Vector4::point(0., 2.0_f32.sqrt() / 2.0, -2.0_f32.sqrt() / 2.0);
        let normal = sphere.normal(point);
        assert_relative_eq!(
            normal,
            Vector4::vector(0., 0.97014, -0.24254),
            epsilon = EPSILON
        );
    }
}
