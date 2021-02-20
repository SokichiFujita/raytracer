use crate::shape::Shape;
use na::Vector4;

#[derive(Clone, Debug)]
pub struct Computations {
    pub t: f32,
    pub shape: Shape,
    pub point: Vector4<f32>,
    pub under_point: Vector4<f32>,
    pub over_point: Vector4<f32>,
    pub eyev: Vector4<f32>,
    pub normalv: Vector4<f32>,
    pub inside: bool,
    pub reflectv: Vector4<f32>,
    pub n1: f32,
    pub n2: f32,
}

impl Computations {
    pub fn schlick(&self) -> f32 {
        let mut cos = self.eyev.dot(&self.normalv);
        if self.n1 > self.n2 {
            let n = self.n1 / self.n2;
            let sin2_t = n.powi(2) * (1.0 - cos.powi(2));
            if sin2_t > 1.0 {
                return 1.0;
            }
            let cos_t = (1.0 - sin2_t).sqrt();
            cos = cos_t;
        }
        let r0 = ((self.n1 - self.n2) / (self.n1 + self.n2)).powi(2);
        r0 + ((1.0 - r0) * (1.0 - cos).powi(5))
    }
}

#[cfg(test)]
mod tests {}
