use crate::tuple::TupleOperation;
use nalgebra::{Matrix4, Vector4};

#[derive(Debug, Clone)]
pub struct Ray {
    pub origin: Vector4<f32>,
    pub direction: Vector4<f32>,
}

impl Ray {
    pub fn new(origin: Vector4<f32>, direction: Vector4<f32>) -> Self {
        Ray { origin, direction }
    }

    pub fn from_tuple(origin: (f32, f32, f32), direction: (f32, f32, f32)) -> Self {
        Ray {
            origin: Vector4::point(origin.0, origin.1, origin.2),
            direction: Vector4::vector(direction.0, direction.1, direction.2),
        }
    }

    pub fn position(&self, distance: f32) -> Vector4<f32> {
        self.origin + self.direction * distance
    }

    pub fn transform(&self, transform_matrix: Matrix4<f32>) -> Ray {
        Self::new(
            transform_matrix * self.origin,
            transform_matrix * self.direction,
        )
    }
    pub fn inv_transform(&self, transform_matrix: Matrix4<f32>) -> Ray {
        Self::new(
            transform_matrix.try_inverse().unwrap() * self.origin,
            transform_matrix.try_inverse().unwrap() * self.direction,
        )
    }
    pub fn reflect(&self, normal_vector: Vector4<f32>) -> Vector4<f32> {
        -(normal_vector * 2.0 * Vector4::dot(&self.direction, &normal_vector) - self.direction)
    }
}

#[cfg(test)]
mod tests {
    use crate::matrix::CGMatrix;
    use crate::{shape::Shape, sphere::Sphere};

    use super::*;
    use nalgebra::Vector4;

    //Ch5 P.58
    #[test]
    fn create_ray() {
        // creating and querying a ray
        let origin = Vector4::point(1.0, 2.0, 3.0);
        let direction = Vector4::vector(4.0, 5.0, 6.0);
        let ray = Ray::new(origin, direction);
        assert_eq!(ray.origin, origin);
        assert_eq!(ray.direction, direction);
    }

    //Ch5 P.58
    #[test]
    fn computing_point() {
        // creating and querying a ray
        let ray = Ray::from_tuple((2.0, 3.0, 4.0), (1.0, 0.0, 0.0));
        assert_eq!(ray.position(0.0), Vector4::point(2.0, 3.0, 4.0));
        assert_eq!(ray.position(1.0), Vector4::point(3.0, 3.0, 4.0));
        assert_eq!(ray.position(-1.0), Vector4::point(1.0, 3.0, 4.0));
        assert_eq!(ray.position(2.5), Vector4::point(4.5, 3.0, 4.0));
    }

    //P.69
    #[test]
    fn translating_ray() {
        // creating and querying a ray
        let ray = Ray::from_tuple((1.0, 2.0, 3.0), (0.0, 1.0, 0.0));
        let translation = Matrix4::translation(3.0, 4.0, 5.0);
        let ray2 = ray.transform(translation);
        assert_eq!(ray2.origin, Vector4::point(4.0, 6.0, 8.0));
        assert_eq!(ray2.direction, Vector4::vector(0.0, 1.0, 0.0));
    }

    //P.69
    #[test]
    fn scaling_ray() {
        // creating and querying a ray
        let ray = Ray::from_tuple((1.0, 2.0, 3.0), (0.0, 1.0, 0.0));
        let scaling = Matrix4::scaling(2.0, 3.0, 4.0);
        let ray2 = ray.transform(scaling);
        assert_eq!(ray2.origin, Vector4::point(2.0, 6.0, 12.0));
        assert_eq!(ray2.direction, Vector4::vector(0.0, 3.0, 0.0));
    }

    //P.69
    #[test]
    fn intersecting_scaled_sphere_ray() {
        // creating and querying a ray
        let ray = Ray::from_tuple((0.0, 0.0, -5.0), (0.0, 0.0, 1.0));
        let s = Shape::Sphere(Sphere::new(
            Some(Matrix4::<f32>::scaling(2., 2., 2.)),
            None,
            None,
            None,
            None,
        ));
        let xs = s.intersect(&ray);
        assert_eq!(xs.len(), 2);
        assert_eq!(xs[0].t, 3.0);
        assert_eq!(xs[1].t, 7.0);
    }
}
