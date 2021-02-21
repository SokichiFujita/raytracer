use crate::{
    cone::Cone, cube::Cube, cylinder::Cylinder, epsilon::EPSILON, intersection::Intersection,
    material::Material, matrix::CGMatrix, plane::Plane, ray::Ray, sphere::Sphere,
    tuple::TupleOperation,
};
use na::{Matrix4, Vector4};
use std::mem::swap;
use ulid::Ulid;

#[derive(Clone, Debug)]
pub enum Shape {
    Sphere(Sphere),
    Plane(Plane),
    Cube(Cube),
    Cylinder(Cylinder),
    Cone(Cone),
}

impl PartialEq for Shape {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Shape::Sphere(x), Shape::Sphere(y)) => x.id == y.id,
            (Shape::Plane(x), Shape::Plane(y)) => x.id == y.id,
            (Shape::Cube(x), Shape::Cube(y)) => x.id == y.id,
            (Shape::Cylinder(x), Shape::Cylinder(y)) => x.id == y.id,
            (Shape::Cone(x), Shape::Cone(y)) => x.id == y.id,
            _ => false,
        }
    }
}

#[allow(dead_code)]
impl Shape {
    pub fn new_sphere() -> Self {
        Shape::Sphere(Sphere::new(None, None, None, None, None))
    }
    pub fn new_sphere_default() -> Self {
        Shape::Sphere(Sphere::new_default())
    }
    pub fn new_plane() -> Self {
        Shape::Plane(Plane::new(None, None, None, None))
    }
    pub fn new_plane_default() -> Self {
        Shape::Plane(Plane::new_default())
    }
    pub fn new_world_defaults() -> Vec<Self> {
        vec![
            Shape::Sphere(Sphere::new(
                None,
                Some(Material::new_default_world()),
                None,
                None,
                None,
            )),
            Shape::Sphere(Sphere::new(
                Some(Matrix4::scaling(0.5, 0.5, 0.5)),
                None,
                None,
                None,
                None,
            )),
        ]
    }

    pub fn generate_id(prefix: Option<&str>) -> String {
        let ulid = Ulid::new().to_string();
        if prefix.is_some() {
            return format!("{}-{}", prefix.unwrap(), ulid).to_string();
        }
        ulid
    }

    pub fn set_material(&mut self, material: Material) -> Self {
        match self {
            Shape::Plane(x) => {
                x.material = material.clone();
                Shape::Plane(x.clone())
            }
            Shape::Sphere(x) => {
                x.material = material.clone();
                Shape::Sphere(x.clone())
            }
            Shape::Cube(x) => {
                x.material = material.clone();
                Shape::Cube(x.clone())
            }
            Shape::Cylinder(x) => {
                x.material = material.clone();
                Shape::Cylinder(x.clone())
            }
            Shape::Cone(x) => {
                x.material = material.clone();
                Shape::Cone(x.clone())
            }
        }
    }

    pub fn normal(&self, point: Vector4<f32>) -> Vector4<f32> {
        match self {
            Shape::Plane(x) => {
                x.transformation.try_inverse().unwrap() * Vector4::vector(0.0, 1.0, 0.0)
            }
            Shape::Sphere(x) => {
                let point_on_shape = x.transformation.try_inverse().unwrap() * point;
                let normal_at_local = point_on_shape - Vector4::vector(0.0, 0.0, 0.0);
                let normal_at_world =
                    x.transformation.try_inverse().unwrap().transpose() * normal_at_local;

                Vector4::vector(normal_at_world.x, normal_at_world.y, normal_at_world.z).normalize()
            }
            Shape::Cube(x) => {
                let point_on_shape = x.transformation.try_inverse().unwrap() * point;
                let maxc = point_on_shape
                    .x
                    .abs()
                    .max(point_on_shape.y.abs())
                    .max(point_on_shape.z.abs());
                if maxc == point.x.abs() {
                    Vector4::vector(point_on_shape.x, 0., 0.)
                } else if maxc == point_on_shape.y.abs() {
                    Vector4::vector(0., point_on_shape.y, 0.)
                } else {
                    Vector4::vector(0., 0., point_on_shape.z)
                }
            }
            Shape::Cylinder(x) => {
                let point_on_shape = x.transformation.try_inverse().unwrap() * point;
                let dist = point_on_shape.x.powi(2) + point_on_shape.z.powi(2);
                if dist < 1.0 && point_on_shape.y >= (x.max - EPSILON) {
                    Vector4::vector(0.0, 1.0, 0.0)
                } else if dist < 1.0 && point_on_shape.y <= (x.min + EPSILON) {
                    Vector4::vector(0.0, -1.0, 0.0)
                } else {
                    Vector4::vector(point_on_shape.x, 0.0, point_on_shape.z)
                }
            }
            Shape::Cone(x) => {
                let point_on_shape = x.transformation.try_inverse().unwrap() * point;
                let dist = point_on_shape.x.powi(2) + point_on_shape.z.powi(2);
                if dist < x.max.powi(2) && point_on_shape.y >= x.max - EPSILON {
                    Vector4::vector(0.0, 1.0, 0.0)
                } else if dist < x.min.powi(2) && point_on_shape.y <= x.min + EPSILON {
                    Vector4::vector(0.0, -1.0, 0.0)
                } else {
                    Vector4::vector(
                        point_on_shape.x,
                        if point_on_shape.y > 0.0 {
                            -dist.sqrt()
                        } else {
                            dist.sqrt()
                        },
                        point_on_shape.z,
                    )
                }
            }
        }
    }

    pub fn intersect(&self, ray: &Ray) -> Vec<Intersection> {
        match self {
            Shape::Plane(x) => {
                let inv_ray = ray.inv_transform(x.transformation);

                if inv_ray.direction.y.abs() < EPSILON {
                    return vec![];
                }
                let t = -inv_ray.origin.y / inv_ray.direction.y;
                let intersection = Intersection::new(t, &self);
                vec![intersection]
            }
            Shape::Sphere(x) => {
                let inv_ray = ray.inv_transform(x.transformation);

                let sphere_to_ray = inv_ray.origin - Vector4::point(0.0, 0.0, 0.0);
                let a = inv_ray.direction.dot(&inv_ray.direction);
                let b = 2.0 * inv_ray.direction.dot(&sphere_to_ray);
                let c = &sphere_to_ray.dot(&sphere_to_ray) - 1.0;
                let discriminant = b.powi(2) - 4.0 * a * c;
                let d = discriminant.sqrt();
                if discriminant < 0.0 {
                    return vec![];
                }
                let t1 = (-b - d) / (2.0 * a);
                let t2 = (-b + d) / (2.0 * a);
                let i1 = Intersection::new(t1, &self);
                let i2 = Intersection::new(t2, &self);
                vec![i1, i2]
            }
            Shape::Cube(x) => {
                let inv_ray = ray.inv_transform(x.transformation);
                let (xtmin, xtmax) = Self::check_cube_axis(inv_ray.origin.x, inv_ray.direction.x);
                let (ytmin, ytmax) = Self::check_cube_axis(inv_ray.origin.y, inv_ray.direction.y);
                let (ztmin, ztmax) = Self::check_cube_axis(inv_ray.origin.z, inv_ray.direction.z);

                let tmin = xtmin.max(ytmin).max(ztmin);
                let tmax = xtmax.min(ytmax).min(ztmax);
                if tmin > tmax {
                    return vec![];
                }
                vec![
                    Intersection::new(tmin, &self),
                    Intersection::new(tmax, &self),
                ]
            }
            Shape::Cylinder(x) => {
                let inv_ray = ray.inv_transform(x.transformation);

                let mut intersections: Vec<Intersection> = vec![];
                let a = inv_ray.direction.x.powi(2) + inv_ray.direction.z.powi(2);
                if approx::relative_eq!(a, 0., epsilon = EPSILON) {
                    let xs = Self::intersect_caps(&self, &inv_ray);
                    return xs;
                } else {
                    let b = 2. * inv_ray.origin.x * inv_ray.direction.x
                        + 2. * inv_ray.origin.z * inv_ray.direction.z;
                    let c = inv_ray.origin.x.powi(2) + inv_ray.origin.z.powi(2) - 1.;
                    let d = b.powi(2) - 4. * a * c;
                    if d < 0. {
                        return vec![];
                    }

                    let mut t0 = (-b - d.sqrt()) / (2. * a);
                    let mut t1 = (-b + d.sqrt()) / (2. * a);
                    if t0 > t1 {
                        swap(&mut t0, &mut t1);
                    }

                    let y0 = inv_ray.origin.y + t0 * inv_ray.direction.y;
                    if x.min < y0 && y0 < x.max {
                        intersections.push(Intersection::new(t0, &self));
                    }

                    let y1 = inv_ray.origin.y + t1 * inv_ray.direction.y;
                    if x.min < y1 && y1 < x.max {
                        intersections.push(Intersection::new(t1, &self));
                    }
                    let xs = [
                        intersections.as_slice(),
                        Self::intersect_caps(&self, &inv_ray).as_slice(),
                    ]
                    .concat();
                    xs
                }
            }
            Shape::Cone(x) => {
                let inv_ray = ray.inv_transform(x.transformation);

                let mut intersections: Vec<Intersection> = vec![];

                let a = inv_ray.direction.x.powi(2) - inv_ray.direction.y.powi(2)
                    + inv_ray.direction.z.powi(2);
                let b = 2.0
                    * (inv_ray.origin.x * inv_ray.direction.x
                        - inv_ray.origin.y * inv_ray.direction.y
                        + inv_ray.origin.z * inv_ray.direction.z);
                let c =
                    inv_ray.origin.x.powi(2) - inv_ray.origin.y.powi(2) + inv_ray.origin.z.powi(2);

                if approx::relative_eq!(a, 0., epsilon = EPSILON)
                    && !approx::relative_eq!(b, 0., epsilon = EPSILON)
                {
                    intersections.push(Intersection::new(c / (-2.0 * b), &self));
                } else {
                    let d = ((b.powi(2) - 4. * a * c) * 100.).round() / 100.; // for avoiding negative small number near zero
                    if d < 0. {
                        return vec![];
                    }

                    let mut t0 = (-b - d.sqrt()) / (2. * a);
                    let mut t1 = (-b + d.sqrt()) / (2. * a);
                    if t0 > t1 {
                        swap(&mut t0, &mut t1);
                    }

                    let y0 = inv_ray.origin.y + t0 * inv_ray.direction.y;
                    if x.min < y0 && y0 < x.max {
                        intersections.push(Intersection::new(t0, &self));
                    }

                    let y1 = inv_ray.origin.y + t1 * inv_ray.direction.y;
                    if x.min < y1 && y1 < x.max {
                        intersections.push(Intersection::new(t1, &self));
                    }
                }
                let xs = [
                    intersections.as_slice(),
                    Self::intersect_caps(&self, &inv_ray).as_slice(),
                ]
                .concat();
                xs
            }
        }
    }
    //Helper for cube
    fn check_cube_axis(origin: f32, direction: f32) -> (f32, f32) {
        let t0: f32;
        let t1: f32;
        if direction >= 0. {
            t0 = (-1. - origin) / direction;
            t1 = (1. - origin) / direction;
        } else {
            t1 = (-1. - origin) / direction;
            t0 = (1. - origin) / direction;
        }
        (t0, t1)
    }

    fn check_cylinder_cap(ray: &Ray, t: f32) -> bool {
        let x = ray.origin.x + t * ray.direction.x;
        let z = ray.origin.z + t * ray.direction.z;
        x.powi(2) + z.powi(2) <= 1.0 + EPSILON
    }

    fn check_cone_cap(ray: &Ray, t: f32, r: f32) -> bool {
        let x = ray.origin.x + t * ray.direction.x;
        let z = ray.origin.z + t * ray.direction.z;
        x.powi(2) + z.powi(2) <= (r).powi(2) + EPSILON
    }

    fn intersect_caps(&self, ray: &Ray) -> Vec<Intersection> {
        let mut intersections = vec![];
        match self {
            Shape::Cylinder(x) => {
                if !x.closed || approx::relative_eq!(ray.direction.y, 0., epsilon = EPSILON) {
                    return vec![];
                }
                let t_lower = (x.min - ray.origin.y) / ray.direction.y;
                if Shape::check_cylinder_cap(&ray, t_lower) {
                    intersections.push(Intersection::new(t_lower, &self));
                }

                let t_upper = (x.max - ray.origin.y) / ray.direction.y;
                if Shape::check_cylinder_cap(&ray, t_upper) {
                    intersections.push(Intersection::new(t_upper, &self));
                }
            }
            Shape::Cone(x) => {
                if !x.closed || approx::relative_eq!(ray.direction.y, 0., epsilon = EPSILON) {
                    return vec![];
                }
                let t_lower = (x.min - ray.origin.y) / ray.direction.y;
                if Shape::check_cone_cap(&ray, t_lower, x.min) {
                    intersections.push(Intersection::new(t_lower, &self));
                }

                let t_upper = (x.max - ray.origin.y) / ray.direction.y;
                if Shape::check_cone_cap(&ray, t_upper, x.max) {
                    intersections.push(Intersection::new(t_upper, &self));
                }
            }
            _ => return vec![],
        }
        intersections
    }

    pub fn material(&self) -> &Material {
        match self {
            Shape::Plane(x) => &x.material,
            Shape::Sphere(x) => &x.material,
            Shape::Cube(x) => &x.material,
            Shape::Cylinder(x) => &x.material,
            Shape::Cone(x) => &x.material,
        }
    }
    pub fn plane(&self) -> Option<&Plane> {
        match self {
            Shape::Plane(x) => Some(x),
            _ => None,
        }
    }
    pub fn sphere(&self) -> Option<&Sphere> {
        match self {
            Shape::Sphere(x) => Some(x),
            _ => None,
        }
    }
    pub fn cube(&self) -> Option<&Cube> {
        match self {
            Shape::Cube(x) => Some(x),
            _ => None,
        }
    }
    pub fn cylinder(&self) -> Option<&Cylinder> {
        match self {
            Shape::Cylinder(x) => Some(x),
            _ => None,
        }
    }
    pub fn cone(&self) -> Option<&Cone> {
        match self {
            Shape::Cone(x) => Some(x),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {}
