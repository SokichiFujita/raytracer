use crate::color::Color;
use crate::{shape::Shape, tuple::TupleOperation};
use na::{Matrix4, Vector4};

#[derive(Copy, Clone, PartialEq, Debug)]
pub enum Pattern {
    Stripe(Stripe),
    Gradient(Gradient),
    Ring(Ring),
    Checker(Checker),
    Test(Test),
}

#[derive(Clone, PartialEq, Debug, Copy)]
pub struct Stripe {
    transformation: Matrix4<f32>,
    c1: Color,
    c2: Color,
}
#[derive(Clone, PartialEq, Debug, Copy)]
pub struct Gradient {
    transformation: Matrix4<f32>,
    c1: Color,
    c2: Color,
}
#[derive(Clone, PartialEq, Debug, Copy)]
pub struct Ring {
    transformation: Matrix4<f32>,
    c1: Color,
    c2: Color,
}

#[derive(Clone, PartialEq, Debug, Copy)]
pub struct Checker {
    transformation: Matrix4<f32>,
    c1: Color,
    c2: Color,
}

#[derive(Clone, PartialEq, Debug, Copy)]
pub struct Test {
    transformation: Matrix4<f32>,
}

impl Test {
    pub fn new() -> Test {
        Test {
            transformation: Matrix4::<f32>::identity(),
        }
    }
}

impl Stripe {
    pub fn new(c1: (f32, f32, f32), c2: (f32, f32, f32), transformation: Matrix4<f32>) -> Stripe {
        Stripe {
            transformation,
            c1: Color::from_tuple(c1),
            c2: Color::from_tuple(c2),
        }
    }
}

impl Gradient {
    pub fn new(c1: (f32, f32, f32), c2: (f32, f32, f32), transformation: Matrix4<f32>) -> Gradient {
        Gradient {
            transformation,
            c1: Color::from_tuple(c1),
            c2: Color::from_tuple(c2),
        }
    }
}

impl Ring {
    pub fn new(c1: (f32, f32, f32), c2: (f32, f32, f32), transformation: Matrix4<f32>) -> Ring {
        Ring {
            transformation,
            c1: Color::from_tuple(c1),
            c2: Color::from_tuple(c2),
        }
    }
}

impl Checker {
    pub fn new(c1: (f32, f32, f32), c2: (f32, f32, f32), transformation: Matrix4<f32>) -> Checker {
        Checker {
            transformation,
            c1: Color::from_tuple(c1),
            c2: Color::from_tuple(c2),
        }
    }
}

impl Pattern {
    pub fn pattern_at(&self, _shape: &Shape, point: Vector4<f32>) -> Color {
        match self {
            Pattern::Stripe(pattern) => {
                let p = point.x.floor();
                let color = if p % 2.0 == 0.0 {
                    pattern.c1
                } else {
                    pattern.c2
                };
                color
            }
            Pattern::Gradient(pattern) => (pattern.c1.vector()
                + (point.x - point.x.floor()) * (pattern.c2 - pattern.c1).vector())
            .to_color(),
            Pattern::Ring(pattern) => {
                let d = ((point.x.powi(2) + point.z.powi(2)).sqrt() % 2.0).floor();
                let color = if d == 0. { pattern.c1 } else { pattern.c2 };
                color
            }
            Pattern::Checker(pattern) => {
                let d = (point.x.floor() + point.y.floor() + point.z.floor()) % 2.0;
                let color = if d == 0. { pattern.c1 } else { pattern.c2 };
                color
            }
            Pattern::Test(_pattern) => Color::rgb(point.x, point.y, point.z),
        }
    }
}

#[cfg(test)]
mod tests {}
