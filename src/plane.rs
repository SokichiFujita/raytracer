use na::{Matrix4, Vector4};
use ulid::Ulid;

use crate::material::Material;
use crate::tuple::TupleOperation;

#[derive(Clone, Debug)]
pub struct Plane {
    pub id: String,
    pub transformation: Matrix4<f32>,
    pub material: Material,
    pub position: Vector4<f32>, //only the y-axis is used. size is infinite.
}

impl PartialEq for Plane {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Plane {
    pub fn new(
        transformation: Option<Matrix4<f32>>,
        matrial: Option<Material>,
        position: Option<Vector4<f32>>,
    ) -> Plane {
        Plane {
            id: Ulid::new().to_string(),
            transformation: match transformation {
                Some(x) => x,
                None => Matrix4::<f32>::identity(),
            },
            material: match matrial {
                Some(x) => x,
                None => Material::new_default(),
            },
            position: match position {
                Some(x) => x,
                None => Vector4::point(0.0, 0.0, 0.0),
            },
        }
    }
    pub fn new_default() -> Plane {
        Plane::new(None, None, None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ray::Ray, shape::Shape};

    //P.123
    #[test]
    fn intersecting_ray_parallel_to_plane() {
        let p = Shape::Plane(Plane::new_default());
        let ray = Ray::from_tuple((0., 10., 0.), (0., 0., 1.));
        let xs = p.intersect(&ray);
        assert_eq!(xs.len(), 0);
    }

    //P.123
    #[test]
    fn intersecting_coplanar_ray() {
        let p = Shape::Plane(Plane::new_default());
        let ray = Ray::from_tuple((0., 0., 0.), (0., 0., 1.));
        let xs = p.intersect(&ray);
        assert_eq!(xs.len(), 0);
    }

    //P.123
    #[test]
    fn intersecting_from_above() {
        let p = Shape::Plane(Plane::new_default());
        let ray = Ray::from_tuple((0., 1., 0.), (0., -1., 0.));
        let xs = p.intersect(&ray);
        assert_eq!(xs.len(), 1);
        assert_eq!(xs.get(0).unwrap().shape, &p);
    }

    //P.123
    #[test]
    fn intersecting_from_below() {
        let p = Shape::Plane(Plane::new_default());
        let ray = Ray::from_tuple((0., -1., 0.), (0., -1., 0.));
        let xs = p.intersect(&ray);
        assert_eq!(xs.len(), 1);
        assert_eq!(xs.get(0).unwrap().shape, &p);
    }
}
