use crate::epsilon::EPSILON;
use approx::RelativeEq;
use core::ops;
use na::Vector4;

#[derive(Clone, Debug, Copy)]
pub struct Color {
    pub r: f32,
    pub g: f32,
    pub b: f32,
}
pub struct ColorU8 {
    pub r: u8,
    pub g: u8,
    pub b: u8,
}

impl ops::Add<Color> for Color {
    type Output = Color;
    fn add(self, other: Color) -> Color {
        Color {
            r: self.r + other.r,
            g: self.g + other.g,
            b: self.b + other.b,
        }
    }
}

impl ops::Sub<Color> for Color {
    type Output = Color;
    fn sub(self, other: Color) -> Color {
        Color {
            r: self.r - other.r,
            g: self.g - other.g,
            b: self.b - other.b,
        }
    }
}

impl ops::Neg for Color {
    type Output = Color;
    fn neg(self) -> Color {
        Color {
            r: -self.r,
            g: -self.g,
            b: -self.b,
        }
    }
}

impl ops::Mul<f32> for Color {
    type Output = Color;
    fn mul(self, k: f32) -> Color {
        Color {
            r: self.r * k,
            g: self.g * k,
            b: self.b * k,
        }
    }
}

impl PartialEq for Color {
    fn eq(&self, other: &Self) -> bool {
        let epsilon = EPSILON;
        f32::relative_eq(&self.r, &other.r, epsilon, epsilon)
            && f32::relative_eq(&self.g, &other.g, epsilon, epsilon)
            && f32::relative_eq(&self.b, &other.b, epsilon, epsilon)
    }
}

impl Color {
    pub fn from_tuple(rgb: (f32, f32, f32)) -> Color {
        Color {
            r: rgb.0,
            g: rgb.1,
            b: rgb.2,
        }
    }

    pub fn rgb(r: f32, g: f32, b: f32) -> Color {
        Color { r, g, b }
    }
    pub fn vector(&self) -> Vector4<f32> {
        Vector4::new(self.r, self.g, self.b, 0.0)
    }
    pub fn from_vector(vector: Vector4<f32>) -> Color {
        Color::rgb(vector.x, vector.y, vector.z)
    }

    pub fn hadamard_product(&self, other: &Color) -> Color {
        Color {
            r: self.r * other.r,
            g: self.g * other.g,
            b: self.b * other.b,
        }
    }
    pub fn to_u8(&self) -> ColorU8 {
        let r = Color::regularize((self.r * 255.0).round());
        let g = Color::regularize((self.g * 255.0).round());
        let b = Color::regularize((self.b * 255.0).round());
        ColorU8 { r, g, b }
    }
    fn regularize(c: f32) -> u8 {
        if c < 0.0 {
            return 0 as u8;
        }
        if c > 255.0 {
            return 255 as u8;
        }
        c as u8
    }
}
