use crate::{material::Material, pointlight, shape::Shape, tuple::TupleOperation};
use approx::RelativeEq;
use core::f32;
use na::{Matrix4, Vector4};
use std::f32::INFINITY;

#[derive(Clone, Debug)]
pub struct Bound {
    pub id: String,
    pub min: Vector4<f32>,
    pub max: Vector4<f32>,
}

impl PartialEq for Bound {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Bound {
    pub fn new(min: Option<Vector4<f32>>, max: Option<Vector4<f32>>) -> Bound {
        Bound {
            id: Shape::generate_id(Some("bound")),
            min: match min {
                Some(x) => x,
                None => Vector4::point(INFINITY, INFINITY, INFINITY),
            },
            max: match max {
                Some(x) => x,
                None => Vector4::point(-INFINITY, -INFINITY, -INFINITY),
            },
        }
    }
    pub fn new_default() -> Bound {
        Bound::new(None, None)
    }
    pub fn new_from_tuple(min: (f32, f32, f32), max: (f32, f32, f32)) -> Bound {
        Bound::new(
            Some(Vector4::point(min.0, min.1, min.2)),
            Some(Vector4::point(max.0, max.1, max.2)),
        )
    }

    pub fn update_point_of_bound(&mut self, point: Vector4<f32>) {
        if point.x < self.min.x {
            self.min.x = point.x
        };
        if point.y < self.min.y {
            self.min.y = point.y
        };
        if point.z < self.min.z {
            self.min.z = point.z
        };
        if point.x > self.max.x {
            self.max.x = point.x
        };
        if point.y > self.max.y {
            self.max.y = point.y
        };
        if point.z > self.max.z {
            self.max.z = point.z
        };
    }

    pub fn update_bound(&mut self, other: Self) {
        self.update_point_of_bound(other.min);
        self.update_point_of_bound(other.max);
    }

    pub fn point_is_contained(&self, point: Vector4<f32>) -> bool {
        let x_is_contained = approx::relative_eq!(self.min.x, point.x)
            || approx::relative_eq!(self.max.x, point.x)
            || (self.min.x < point.x && point.x < self.max.x);
        let y_is_contained = approx::relative_eq!(self.min.y, point.y)
            || approx::relative_eq!(self.max.y, point.y)
            || (self.min.y < point.y && point.y < self.max.y);
        let z_is_contained = approx::relative_eq!(self.min.z, point.z)
            || approx::relative_eq!(self.max.z, point.z)
            || (self.min.z < point.z && point.z < self.max.z);
        x_is_contained && y_is_contained && z_is_contained
    }

    pub fn bound_is_contained(&self, other: Bound) -> bool {
        self.point_is_contained(other.min) && self.point_is_contained(other.max)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ray::Ray, shape::Shape, tuple::TupleOperation};
    use approx::assert_relative_eq;
    use nalgebra::Vector4;
}
