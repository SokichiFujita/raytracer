use crate::{
    color::Color, intersection::*, pointlight::PointLight, ppm::Ppm, ray::Ray, shape::Shape,
    tuple::TupleOperation,
};
use na::{Vector, Vector4};

#[derive(Clone, Debug)]
pub struct Canvas {
    pub width: usize,
    pub height: usize,
    pub pixels: Vec<Color>,
}

impl Canvas {
    pub fn new(width: usize, height: usize, pixels: Vec<Color>) -> Canvas {
        Canvas {
            width,
            height,
            pixels,
        }
    }

    pub fn to_ppm(&self) -> Ppm {
        Ppm::new(self.width, self.height, &self.pixels)
    }

    pub fn render_single_shape(shape: Shape, width: usize) -> Canvas {
        let canvas_pixels = width;
        let ray_origin_z = -5.0;
        let wall_z = 10.0;
        let wall_size = 7.0;
        let pixel_size = wall_size / (width as f32);
        let half = wall_size / 2.0;
        let ray = Ray::from_tuple((0.0, 0.0, ray_origin_z), (0.0, 0.0, wall_z));
        let light = PointLight::from_tuple((-10., -10., -10.), (1., 1., 1.));

        let mut pixels: Vec<Color> = vec![];

        for y in (0..canvas_pixels).rev() {
            let world_y = half - pixel_size * (y as f32);
            for x in 0..canvas_pixels {
                let world_x = -half + pixel_size * (x as f32);
                let position = Vector4::point(world_x, world_y, wall_z);
                let direction = Vector::normalize(&(position - ray.origin));
                let r = Ray::new(ray.origin, direction);
                let intersections = shape.intersect(&r);
                let hit = intersections.hit();
                if hit.is_some() {
                    let point = r.position(hit.unwrap().t);
                    let normalv = hit.unwrap().shape.normal(point, None);
                    let eyev = -r.direction;
                    let shape = hit.unwrap().shape;
                    match shape {
                        Shape::Sphere(x) => {
                            let color = x.material.lighting(
                                &light,
                                point,
                                eyev,
                                normalv,
                                false,
                                Some(&shape),
                            );
                            pixels.push(color);
                        }
                        _ => {}
                    }
                } else {
                    let black = Color::rgb(0.0, 0.0, 0.0);
                    pixels.push(black);
                }
            }
        }
        Canvas {
            width: width,
            height: width,
            pixels,
        }
    }
}

#[cfg(test)]
mod tests {}
