use crate::{
    cone::Cone, cube::Cube, cylinder::Cylinder, epsilon::EPSILON, intersection::Intersection,
    material::Material, matrix::CGMatrix, plane::Plane, ray::Ray, shape::Shape, sphere::Sphere,
    tuple::TupleOperation,
};
use na::{Matrix4, Vector4};
use std::mem::swap;
use ulid::Ulid;

#[derive(Clone, Debug)]
pub struct Group {
    pub id: String,
    pub transformation: Matrix4<f32>,
    pub material: Material,
    pub children: Vec<Shape>,
    pub parent: Option<String>,
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
}

#[cfg(test)]
mod tests {}
