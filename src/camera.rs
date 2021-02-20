use na::{Matrix4, Vector4};

use crate::canvas::Canvas;
use crate::color::Color;
use crate::ray::Ray;
use crate::world::World;

use crate::tuple::TupleOperation;

#[derive(Clone, Debug)]
pub struct Camera {
    horizontal_size: f32,
    vertical_size: f32,
    field_of_view: f32,
    transform: Matrix4<f32>,
    half_width: f32,
    half_height: f32,
    pixel_size: f32,
}

#[allow(dead_code)]
impl Camera {
    pub fn new(
        horizontal_size: f32,
        vertical_size: f32,
        field_of_view: f32,
        transform: Matrix4<f32>,
    ) -> Camera {
        let half_view = (field_of_view / 2.0).tan();
        let aspect_ratio = horizontal_size / vertical_size;
        let half_width = if aspect_ratio >= 1.0 {
            half_view
        } else {
            half_view * aspect_ratio
        };
        let half_height = if aspect_ratio >= 1.0 {
            half_view / aspect_ratio
        } else {
            half_view
        };
        let pixel_size = (half_width * 2.0) / horizontal_size;
        Camera {
            horizontal_size,
            vertical_size,
            field_of_view,
            transform,
            half_width,
            half_height,
            pixel_size,
        }
    }

    pub fn new_default(horizontal_size: f32, vertical_size: f32, field_of_view: f32) -> Camera {
        Camera::new(
            horizontal_size,
            vertical_size,
            field_of_view,
            Matrix4::identity(),
        )
    }

    pub fn ray_for_pixel(&self, px: f32, py: f32) -> Ray {
        let x_offset = (px + 0.5) * self.pixel_size;
        let y_offset = (py + 0.5) * self.pixel_size;
        let world_x = self.half_width - x_offset;
        let world_y = self.half_height - y_offset;
        let inv = self.transform.try_inverse().unwrap();
        let pixel = inv * Vector4::point(world_x, world_y, -1.0);
        let origin = inv * Vector4::point(0.0, 0.0, 0.0);
        let direction = (pixel - origin).normalize();
        Ray::new(origin, direction)
    }

    pub fn render(&self, world: &World) -> Canvas {
        let mut pixels = vec![] as Vec<Color>;
        for y in 0..(self.vertical_size as usize) {
            for x in 0..(self.horizontal_size as usize) {
                let ray = self.ray_for_pixel(x as f32, y as f32);
                let color = world.color_at(&ray, 5);
                pixels.push(color);
            }
        }
        Canvas::new(
            self.horizontal_size as usize,
            self.vertical_size as usize,
            pixels,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use std::f32::consts::PI;

    use crate::epsilon::EPSILON;
    use crate::matrix::CGMatrix;
    use crate::world::World;

    // P.95
    #[test]
    fn pixel_size_horizontal_canvas() {
        let hsize = 200.0;
        let vsize = 125.0;
        let field_of_view = PI / 2.0;
        let camera = Camera::new_default(hsize, vsize, field_of_view);
        assert_relative_eq!(camera.pixel_size, 0.01);
    }

    #[test]
    fn pixel_size_vertical_canvas() {
        let hsize = 125.0;
        let vsize = 200.0;
        let field_of_view = PI / 2.0;
        let camera = Camera::new_default(hsize, vsize, field_of_view);
        assert_relative_eq!(camera.pixel_size, 0.01);
    }
    #[test]
    fn constructing_ray_through_center_of_canvas() {
        let camera = Camera::new_default(201.0, 101.0, PI / 2.0);
        let ray = camera.ray_for_pixel(100.0, 50.0);
        assert_relative_eq!(ray.origin, Vector4::new(0.0, 0.0, 0.0, 1.0));
        assert_relative_eq!(ray.direction, Vector4::new(0.0, 0.0, -1.0, 0.0));
    }
    #[test]
    fn constructing_ray_through_corner_of_canvas() {
        let camera = Camera::new_default(201.0, 101.0, PI / 2.0);
        let ray = camera.ray_for_pixel(0.0, 0.0);
        assert_relative_eq!(ray.origin, Vector4::new(0.0, 0.0, 0.0, 1.0));
        assert_relative_eq!(
            ray.direction,
            Vector4::new(0.66519, 0.33259, -0.66851, 0.0),
            epsilon = EPSILON
        );
    }

    #[test]
    fn constructing_ray_camera_is_transformed() {
        let camera = Camera::new(
            201.0,
            101.0,
            PI / 2.0,
            Matrix4::rotation_y(PI / 4.0) * Matrix4::translation(0.0, -2.0, 5.0),
        );
        let ray = camera.ray_for_pixel(100.0, 50.0);
        assert_relative_eq!(ray.origin, Vector4::new(0.0, 2.0, -5.0, 1.0));
        assert_relative_eq!(
            ray.direction,
            Vector4::new(2.0_f32.sqrt() / 2.0, 0.0, -2.0_f32.sqrt() / 2.0, 0.0),
            epsilon = EPSILON
        );
    }

    #[test]
    fn rendering_world_with_camera() {
        let world = World::new_default();
        let transformation = Vector4::point(0.0, 0.0, -5.0).view_transformation(
            &Vector4::point(0.0, 0.0, 0.0),
            &Vector4::vector(0.0, 1.0, 0.0),
        );
        let camera = Camera::new(11.0, 11.0, PI / 2.0, transformation);
        let canvas = camera.render(&world);
        assert_relative_eq!(
            canvas.pixels.get(5 * 11 + 5).unwrap().vector(),
            Color::rgb(0.38066, 0.47583, 0.2855).vector(),
            epsilon = EPSILON
        );
    }
}
