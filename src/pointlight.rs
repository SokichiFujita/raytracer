use na::Vector4;

use crate::color::Color;
use crate::tuple::TupleOperation;

#[derive(Clone, PartialEq, Debug)]
pub struct PointLight {
    pub position: Vector4<f32>,
    pub intensity: Color,
}

#[allow(dead_code)]
impl PointLight {
    pub fn new(position: Vector4<f32>, intensity: Color) -> PointLight {
        PointLight {
            position,
            intensity,
        }
    }
    pub fn new_default() -> PointLight {
        PointLight::from_tuple((-10., 10., -10.), (1., 1., 1.))
    }

    pub fn from_tuple(position: (f32, f32, f32), intensity: (f32, f32, f32)) -> PointLight {
        PointLight {
            position: Vector4::point(position.0, position.1, position.2),
            intensity: Color::rgb(intensity.0, intensity.1, intensity.2),
        }
    }
}

#[cfg(test)]
mod tests {}
