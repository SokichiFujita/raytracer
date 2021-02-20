//use approx::assert_relative_eq;
use na::Matrix4;

pub trait CGMatrix {
    fn rotation_x(x: f32) -> Matrix4<f32>;
    fn rotation_y(x: f32) -> Matrix4<f32>;
    fn rotation_z(x: f32) -> Matrix4<f32>;
    fn translation(x: f32, y: f32, z: f32) -> Matrix4<f32>;
    fn translate(x: f32, y: f32, z: f32) -> Matrix4<f32>;
    fn scaling(x: f32, y: f32, z: f32) -> Matrix4<f32>;
}

//struct CGMatrix;
impl CGMatrix for Matrix4<f32> {
    fn rotation_x(x: f32) -> Matrix4<f32> {
        Matrix4::new(
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            x.cos(),
            -x.sin(),
            0.0,
            0.0,
            x.sin(),
            x.cos(),
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
        )
    }
    fn rotation_y(x: f32) -> Matrix4<f32> {
        Matrix4::new(
            x.cos(),
            0.0,
            x.sin(),
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            -x.sin(),
            0.0,
            x.cos(),
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
        )
    }
    fn rotation_z(x: f32) -> Matrix4<f32> {
        Matrix4::new(
            x.cos(),
            -x.sin(),
            0.0,
            0.0,
            x.sin(),
            x.cos(),
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
        )
    }
    fn scaling(x: f32, y: f32, z: f32) -> Matrix4<f32> {
        Matrix4::new(
            x, 0.0, 0.0, 0.0, 0.0, y, 0.0, 0.0, 0.0, 0.0, z, 0.0, 0.0, 0.0, 0.0, 1.0,
        )
    }
    fn translation(x: f32, y: f32, z: f32) -> Matrix4<f32> {
        Matrix4::new(
            1.0, 0.0, 0.0, x, 0.0, 1.0, 0.0, y, 0.0, 0.0, 1.0, z, 0.0, 0.0, 0.0, 1.0,
        )
    }
    fn translate(x: f32, y: f32, z: f32) -> Matrix4<f32> {
        Matrix4::new(
            1.0, 0.0, 0.0, x, 0.0, 1.0, 0.0, y, 0.0, 0.0, 1.0, z, 0.0, 0.0, 0.0, 1.0,
        )
    }
}
