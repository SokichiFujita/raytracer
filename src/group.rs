use crate::{
    cone::Cone, cube::Cube, cylinder::Cylinder, epsilon::EPSILON, intersection::Intersection,
    material::Material, matrix::CGMatrix, plane::Plane, ray::Ray, shape::Shape, sphere::Sphere,
    tuple::TupleOperation,
};
use na::{Matrix4, Vector4};
use nalgebra::coordinates::M4x4;
use std::mem::swap;

#[derive(Clone, Debug)]
pub struct Group {
    pub id: String,
    pub transformation: Matrix4<f32>,
    pub material: Material,
    pub children: Vec<Shape>,
    pub parent: Option<String>,
}

impl PartialEq for Group {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

#[allow(dead_code)]
impl Group {
    pub fn new(
        transformation: Option<Matrix4<f32>>,
        material: Option<Material>,
        children: Vec<Shape>,
        parent: Option<String>,
    ) -> Group {
        Group {
            id: Shape::generate_id(Some("group")),
            transformation: match transformation {
                Some(x) => x,
                None => Matrix4::<f32>::identity(),
            },
            material: match material {
                Some(x) => x,
                None => Material::new_default(),
            },
            children,
            parent,
        }
    }

    pub fn new_default() -> Group {
        Group::new(None, None, vec![], None)
    }

    pub fn push_to_children(&mut self, shape: &mut Shape) {
        shape.set_parent(self.id.clone());
        self.children.push(shape.clone())
    }

    pub fn intersect(&mut self, ray: &Ray) -> Vec<Intersection> {
        let mut intersections = self
            .children
            .iter()
            .map(|x| {
                let z = x.intersect_group(ray, self.transformation);
                z
            })
            .flatten()
            .collect::<Vec<_>>();
        intersections.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
        intersections
    }
}

#[cfg(test)]
mod tests {
    use std::f32::consts::PI;

    use super::*;
    #[test]
    fn create_new_group() {
        let g = Group::new_default();
        assert_eq!(g.children.len(), 0);
    }

    #[test]
    fn shape_has_parent_attribute() {
        let s = Shape::Plane(Plane::new_default());
        assert_eq!(s.parent(), None);
    }

    #[test]
    fn add_child_to_group() {
        let mut g = Group::new_default();
        let s = Shape::Plane(Plane::new_default());
        g.push_to_children(&mut s.clone());
        assert_eq!(g.children.len(), 1);
        assert_eq!(g.children.get(0).unwrap(), &s);
        assert_eq!(g.children.get(0).unwrap().parent().unwrap(), g.id);
    }

    #[test]
    fn intersectiong_ray_with_nonempty_group() {
        let mut g = Group::new_default();
        let mut s1 = Shape::Sphere(Sphere::new_default());
        let mut s2 = Shape::Sphere(Sphere::new(
            Some(Matrix4::translation(0., 0., -3.)),
            None,
            None,
            None,
            None,
        ));
        let mut s3 = Shape::Sphere(Sphere::new(
            Some(Matrix4::translation(5., 0., 0.)),
            None,
            None,
            None,
            None,
        ));
        g.push_to_children(&mut s1);
        g.push_to_children(&mut s2);
        g.push_to_children(&mut s3);

        let ray = Ray::from_tuple((0., 0., -5.), (0., 0., 1.));
        let xs = g.intersect(&ray);
        println!("{:?}", xs);

        assert_eq!(xs.len(), 4);
        assert_eq!(xs[0].shape, &s2);
        assert_eq!(xs[1].shape, &s2);
        assert_eq!(xs[2].shape, &s1);
        assert_eq!(xs[3].shape, &s1);
    }

    #[test]
    fn intersectiong_transformed_group() {
        let mut g = Group::new(Some(Matrix4::scaling(2., 2., 2.)), None, vec![], None);
        let mut s = Shape::Sphere(Sphere::new(
            Some(Matrix4::translation(5., 0., 0.)),
            None,
            None,
            None,
            None,
        ));
        g.push_to_children(&mut s);

        let ray = Ray::from_tuple((10., 0., -10.), (0., 0., 1.));
        let xs = g.intersect(&ray);
        assert_eq!(xs.len(), 2);
    }

    #[test]
    fn converting_point_from_world_to_object_space() {
        let mut g1 = Group::new(Some(Matrix4::rotation_y(PI / 2.0)), None, vec![], None);
        let mut g2 = Group::new(Some(Matrix4::scaling(2., 2., 2.)), None, vec![], None);
    }
}
