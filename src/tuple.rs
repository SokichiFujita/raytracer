use crate::color::Color;
use crate::matrix::CGMatrix;
use na::{Matrix4, Vector4};

pub trait TupleOperation {
    fn point(x: f32, y: f32, z: f32) -> Vector4<f32>;
    fn vector(x: f32, y: f32, z: f32) -> Vector4<f32>;
    fn is_point(v: Vector4<f32>) -> bool;
    fn is_vector(v: Vector4<f32>) -> bool;
    fn magnitude(&self) -> f32;
    fn reflect(&self, normal: Vector4<f32>) -> Vector4<f32>;
    fn view_transformation(&self, to: &Vector4<f32>, up: &Vector4<f32>) -> Matrix4<f32>;
    fn cross4(&self, other: &Vector4<f32>) -> Vector4<f32>;
    fn to_color(&self) -> Color;
}

impl TupleOperation for Vector4<f32> {
    fn point(x: f32, y: f32, z: f32) -> Vector4<f32> {
        Vector4::new(x, y, z, 1.0)
    }

    fn vector(x: f32, y: f32, z: f32) -> Vector4<f32> {
        Vector4::new(x, y, z, 0.0)
    }

    fn is_point(v: Vector4<f32>) -> bool {
        v[3] == 1.0
    }

    fn is_vector(v: Vector4<f32>) -> bool {
        v[3] == 0.0
    }

    fn magnitude(&self) -> f32 {
        (self.x.powi(2) + self.y.powi(2) + self.z.powi(2) + self[3].powi(2)).sqrt()
    }

    fn reflect(&self, normal: Vector4<f32>) -> Vector4<f32> {
        self - normal * 2.0 * self.dot(&normal)
    }

    fn cross4(&self, other: &Vector4<f32>) -> Vector4<f32> {
        Vector4::vector(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    fn view_transformation(&self, to: &Vector4<f32>, up: &Vector4<f32>) -> Matrix4<f32> {
        let forward = (to - self).normalize();
        let up_normal = up.normalize();
        let left = forward.cross4(&up_normal);
        let true_up = left.cross4(&forward);
        let orientation = Matrix4::<f32>::new(
            left.x, left.y, left.z, 0.0, true_up.x, true_up.y, true_up.z, 0.0, -forward.x,
            -forward.y, -forward.z, 0.0, 0.0, 0.0, 0.0, 1.0,
        );
        orientation * Matrix4::<f32>::translation(-self.x, -self.y, -self.z)
    }

    fn to_color(&self) -> Color {
        Color::rgb(self.x, self.y, self.z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::epsilon::EPSILON;
    use approx::assert_relative_eq;
    //Ch6 P.83
    #[test]
    fn reflecting_a_vector_approaching_at_45() {
        //
        let v = Vector4::vector(1.0, -1.0, 0.0);
        let n = Vector4::vector(0.0, 1.0, 0.0);
        let r = v.reflect(n);
        assert_relative_eq!(r, Vector4::vector(1.0, 1.0, 0.0));
    }

    #[test]
    fn reflecting_a_vector_off_a_slanted_surface() {
        let v = Vector4::vector(1.0, -1.0, 0.0);
        let n = Vector4::vector(0.0, 1.0, 0.0);
        let r = v.reflect(n);
        assert_relative_eq!(r, Vector4::vector(1.0, 1.0, 0.0));
    }

    //ch7 view transformation
    #[test]
    fn transofrmation_matrix_for_default_orientation() {
        let from = Vector4::point(0.0, 0.0, 0.0);
        let to = Vector4::point(0.0, 0.0, -1.0);
        let up = Vector4::vector(0.0, 1.0, 0.0);
        assert_relative_eq!(from.view_transformation(&to, &up), Matrix4::identity());
    }

    //ch7 view transformation
    #[test]
    fn transofrmation_matrix_looking_in_positive_z_direction() {
        let from = Vector4::point(0.0, 0.0, 0.0);
        let to = Vector4::point(0.0, 0.0, 1.0);
        let up = Vector4::vector(0.0, 1.0, 0.0);
        assert_relative_eq!(
            from.view_transformation(&to, &up),
            Matrix4::scaling(-1.0, 1.0, -1.0)
        );
    }

    //ch7 view transformation
    #[test]
    fn transofrmation_moves_world() {
        let from = Vector4::point(0.0, 0.0, 8.0);
        let to = Vector4::point(0.0, 0.0, 0.0);
        let up = Vector4::vector(0.0, 1.0, 0.0);
        assert_relative_eq!(
            from.view_transformation(&to, &up),
            Matrix4::translation(0.0, 0.0, -8.0)
        );
    }

    //ch7 view transformation
    #[test]
    fn arbitrary_view_transformation() {
        let from = Vector4::point(1.0, 3.0, 2.0);
        let to = Vector4::point(4.0, -2.0, 8.0);
        let up = Vector4::vector(1.0, 1.0, 0.0);
        assert_relative_eq!(
            from.view_transformation(&to, &up),
            Matrix4::new(
                -0.50709, 0.50709, 0.67612, -2.36643, 0.76772, 0.60609, 0.12122, -2.82843,
                -0.35857, 0.59761, -0.71714, 0.0, 0.0, 0.0, 0.0, 1.0
            ),
            epsilon = EPSILON
        );
    }
}
