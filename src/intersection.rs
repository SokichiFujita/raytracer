use crate::computation::Computations;
use crate::epsilon::EPSILON;
use crate::ray::Ray;
use crate::shape::*;
use std::cmp::{Ord, Ordering};

use crate::tuple::TupleOperation;

#[derive(Clone, Debug, Copy)]
pub struct Intersection<'a> {
    pub t: f32,
    pub shape: &'a Shape,
    pub u: Option<f32>,
    pub v: Option<f32>,
}

impl PartialEq for Intersection<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.t == other.t && self.shape == other.shape
    }
}

impl PartialOrd for Intersection<'_> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.t.partial_cmp(&other.t)
    }
}

impl Ord for Intersection<'_> {
    fn cmp(&self, other: &Intersection) -> Ordering {
        if self.t < other.t {
            Ordering::Less
        } else if self.t > other.t {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    }
}

impl Eq for Intersection<'_> {}

impl Intersection<'_> {
    pub fn new<'a>(t: f32, shape: &'a Shape) -> Intersection {
        Intersection {
            t,
            shape,
            u: None,
            v: None,
        }
    }
    pub fn new_with_uv<'a>(t: f32, shape: &'a Shape, u: f32, v: f32) -> Intersection {
        Intersection {
            t,
            shape,
            u: Some(u),
            v: Some(v),
        }
    }
    pub fn prepare_computation(
        &self,
        ray: &Ray,
        all_intersections: Option<&Vec<Intersection>>,
    ) -> Computations {
        let point = ray.position(self.t);
        let eyev = -1.0 * ray.direction;
        let hit = match all_intersections {
            None => None,
            Some(x) => Some(x.get(0).unwrap().clone()),
        };
        let normalv_temp = self.shape.normal(point, hit);

        //negating
        let angle_of_normal_and_eye = normalv_temp.dot(&eyev);
        let inside = angle_of_normal_and_eye < 0.0;
        let normalv = if inside {
            -1.0 * normalv_temp
        } else {
            1.0 * normalv_temp
        };

        let offset = EPSILON * normalv;
        let over_point = point + offset;
        let under_point = point - offset;
        let reflectv = ray.direction.reflect(normalv_temp);
        let mut container: Vec<Shape> = vec![];

        let mut n1 = 1.0;
        let mut n2 = 1.0;

        //should be improved
        let mut all: Vec<Intersection> = vec![];
        if all_intersections.is_none() {
            all.push(self.clone());
        } else {
            all = all_intersections.unwrap().clone();
        }

        for x in all.iter() {
            if x == self {
                if container.len() == 0 {
                    n1 = 1.0
                } else {
                    let shape = container.get(container.len() - 1).unwrap();
                    n1 = shape.material().refractive;
                }
            }

            if container.iter().find(|shape| **shape == *x.shape).is_some() {
                container = container
                    .iter()
                    .filter(|shape| **shape != *x.shape)
                    .map(|x| x.clone())
                    .collect::<Vec<_>>();
            } else {
                container.push(x.shape.clone())
            }

            if x == self {
                if container.len() == 0 {
                    n2 = 1.0
                } else {
                    let shape = container.get(container.len() - 1).unwrap();
                    n2 = shape.material().refractive;
                }
            }
        }

        Computations {
            t: self.t,
            shape: self.shape.clone(),
            point,
            under_point,
            over_point,
            eyev,
            normalv,
            inside,
            reflectv,
            n1,
            n2,
        }
    }
}

pub trait Intersections {
    fn hit(&self) -> Option<&Intersection>;
}

impl Intersections for Vec<Intersection<'_>> {
    fn hit(&self) -> Option<&Intersection> {
        let hit = self.iter().filter(|i| i.t >= 0.0).min();
        hit
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        color::Color, material::Material, matrix::CGMatrix, plane::Plane, sphere::Sphere,
        world::World,
    };
    use approx::assert_relative_eq;
    use na::Vector4;
    use nalgebra::Matrix4;

    #[test]
    fn hit_when_all_intersection_positive() {
        let s = Shape::Sphere(Sphere::new_default());
        let i1 = Intersection::new(1.0, &s);
        let i2 = Intersection::new(2.0, &s);
        let xs = vec![i1, i2];
        let i = xs.hit();
        assert_eq!(i, xs.get(0))
    }
    #[test]
    fn hit_when_some_intersection_negative() {
        let s = Shape::Sphere(Sphere::new_default());
        let i1 = Intersection::new(-1.0, &s);
        let i2 = Intersection::new(1.0, &s);
        let xs = vec![i1, i2];
        let i = xs.hit();
        assert_eq!(i, xs.get(1))
    }

    #[test]
    fn hit_when_all_intersection_negative() {
        let s = Shape::Sphere(Sphere::new_default());
        let i1 = Intersection::new(-2.0, &s);
        let i2 = Intersection::new(-1.0, &s);
        let xs = vec![i1, i2];
        let i = xs.hit();
        assert_eq!(i, None)
    }
    #[test]
    fn hit_is_always_lowest_nonnegative_intersection() {
        let s = Shape::Sphere(Sphere::new_default());
        let i1 = Intersection::new(5.0, &s);
        let i2 = Intersection::new(7.0, &s);
        let i3 = Intersection::new(-3.0, &s);
        let i4 = Intersection::new(2.0, &s);
        let xs = vec![i1, i2, i3, i4];
        let i = xs.hit();
        assert_eq!(i, xs.get(3))
    }

    #[test]
    fn precomputing_state_intersection() {
        let ray = Ray::from_tuple((0., 0., -5.), (0., 0., 1.));
        let s = Shape::Sphere(Sphere::new_default());
        let i1 = Intersection::new(4.0, &s);
        let comps = i1.prepare_computation(&ray, None);
        assert_eq!(comps.t, i1.t);
        assert_eq!(comps.shape.clone(), *i1.shape);
        assert_eq!(comps.point, Vector4::point(0., 0., -1.));
        assert_relative_eq!(comps.eyev, Vector4::vector(0., 0., -1.));
        assert_relative_eq!(comps.normalv, Vector4::vector(0., 0., -1.));
    }

    #[test]
    fn hit_of_intersection_occur_outside() {
        let ray = Ray::from_tuple((0., 0., -5.), (0., 0., 1.));
        let s = Shape::Sphere(Sphere::new_default());
        let i1 = Intersection::new(4.0, &s);
        let comps = i1.prepare_computation(&ray, None);
        assert_eq!(comps.inside, false);
    }

    #[test]
    fn hit_of_intersection_occur_inside() {
        let ray = Ray::from_tuple((0., 0., 0.), (0., 0., 1.));
        let s = Shape::Sphere(Sphere::new_default());
        let i1 = Intersection::new(1.0, &s);
        let comps = i1.prepare_computation(&ray, None);
        assert_eq!(comps.inside, true);
        assert_eq!(comps.point, Vector4::point(0., 0., 1.));
        assert_relative_eq!(comps.eyev, Vector4::vector(0., 0., -1.));
        assert_relative_eq!(comps.normalv, Vector4::vector(0., 0., -1.));
    }

    #[test]
    fn hit_should_offset_point() {
        let ray = Ray::from_tuple((0., 0., -5.), (0., 0., 1.));
        let s = Shape::Sphere(Sphere::new(
            Some(Matrix4::<f32>::translation(0., 0., 1.)),
            None,
            None,
            None,
        ));
        let i1 = Intersection::new(5.0, &s);
        let comps = i1.prepare_computation(&ray, None);
        assert_eq!(comps.over_point.z < -EPSILON / 2.0, true);
        assert_eq!(comps.point.z > comps.over_point.z, true);
    }

    #[test]
    fn precomputing_refrection_vector() {
        let ray = Ray::from_tuple(
            (0., 1., -1.),
            (0., -(2_f32.sqrt()) / 2.0, 2_f32.sqrt() / 2.0),
        );
        let s = Shape::Plane(Plane::new_default());
        let i1 = Intersection::new(2_f32.sqrt(), &s);
        let comps = i1.prepare_computation(&ray, None);
        assert_relative_eq!(
            comps.reflectv,
            Vector4::vector(0., 2_f32.sqrt() / 2.0, 2_f32.sqrt() / 2.0)
        );
    }

    #[test]
    fn refrected_color_nonreflective_material() {
        let material = Material::new(
            Some(Color::rgb(0.8, 1.0, 0.6)),
            Some(1.0),
            Some(0.7),
            Some(0.2),
            None,
            None,
            None,
            None,
            None,
        );
        let s1 = Shape::Sphere(Sphere::new(None, Some(material), None, None));
        let s2 = Shape::Sphere(Sphere::new(
            Some(Matrix4::scaling(0.5, 0.5, 0.5)),
            Some(Material::new(
                None,
                Some(1.0),
                None,
                None,
                None,
                None,
                None,
                None,
                None,
            )),
            None,
            None,
        ));
        let ray = Ray::from_tuple((0., 0., 0.), (0., 0., 1.));
        let world = World::new(None, vec![s1.clone(), s2.clone()]);
        let i1 = Intersection::new(1.0, &s2);
        let comps = i1.prepare_computation(&ray, None);
        let color = world.reflected_color(&comps, 1);
        assert_relative_eq!(color.vector(), Color::rgb(0., 0., 0.).vector());
    }

    #[test]
    fn reflected_color_reflective_material() {
        let world = World::new_default();
        let s3 = Shape::Plane(Plane::new(
            Some(Matrix4::<f32>::translate(0., -1., 0.)),
            Some(Material::new(
                None,
                None,
                None,
                None,
                None,
                Some(0.5),
                None,
                None,
                None,
            )),
            None,
        ));
        let ray = Ray::from_tuple(
            (0., 0., -3.),
            (0.0, -(2_f32.sqrt() / 2.0), 2_f32.sqrt() / 2.0),
        );
        let i1 = Intersection::new(2_f32.sqrt(), &s3);
        let comps = i1.prepare_computation(&ray, None);
        let color = world.reflected_color(&comps, 1);
        assert_relative_eq!(
            color.vector(),
            Color::rgb(0.19032, 0.2379, 0.14274).vector(),
            epsilon = EPSILON
        );
    }

    #[test]
    fn reflected_color_maximum_recursive_depth() {
        let mut world = World::new_default();
        let s3 = Shape::Plane(Plane::new(
            Some(Matrix4::<f32>::translation(0., -1., 0.)),
            Some(Material::new(
                None,
                None,
                None,
                None,
                None,
                Some(0.5),
                None,
                None,
                None,
            )),
            None,
        ));
        world.shapes.push(s3.clone());
        let ray = Ray::from_tuple(
            (0., 0., -3.),
            (0.0, -(2_f32.sqrt()) / 2.0, 2_f32.sqrt() / 2.0),
        );
        let i1 = Intersection::new(2_f32.sqrt(), &s3);
        let comps = i1.prepare_computation(&ray, None);
        let color = world.reflected_color(&comps, 0);
        assert_relative_eq!(color.vector(), Color::rgb(0., 0., 0.).vector());
    }

    #[test]
    fn finding_n1_n2() {
        let mut ma = Material::new_glass();
        ma.refractive = 1.5;
        let mut mb = Material::new_glass();
        mb.refractive = 2.0;
        let mut mc = Material::new_glass();
        mc.refractive = 2.5;

        let a = Shape::Sphere(Sphere::new(
            Some(Matrix4::<f32>::scaling(2.0, 2.0, 2.0)),
            Some(ma),
            None,
            None,
        ));
        let b = Shape::Sphere(Sphere::new(
            Some(Matrix4::<f32>::translation(0.0, 0.0, 1.0)),
            Some(mb),
            None,
            None,
        ));
        let c = Shape::Sphere(Sphere::new(
            Some(Matrix4::<f32>::translation(0.0, 0.0, 1.0)),
            Some(mc),
            None,
            None,
        ));

        let ray = Ray::from_tuple((0., 0., -4.), (0., 0., 1.));

        let i1 = Intersection::new(2.0, &a);
        let i2 = Intersection::new(2.75, &b);
        let i3 = Intersection::new(3.25, &c);
        let i4 = Intersection::new(4.75, &b);
        let i5 = Intersection::new(5.25, &c);
        let i6 = Intersection::new(6.0, &a);

        let xs = vec![i1, i2, i3, i4, i5, i6];

        let comps0 = xs[0].prepare_computation(&ray, Some(&xs));
        let comps1 = xs[1].prepare_computation(&ray, Some(&xs));
        let comps2 = xs[2].prepare_computation(&ray, Some(&xs));
        let comps3 = xs[3].prepare_computation(&ray, Some(&xs));
        let comps4 = xs[4].prepare_computation(&ray, Some(&xs));
        let comps5 = xs[5].prepare_computation(&ray, Some(&xs));

        assert_eq!((comps0.n1, comps0.n2), (1.0, 1.5));
        assert_eq!((comps1.n1, comps1.n2), (1.5, 2.0));
        assert_eq!((comps2.n1, comps2.n2), (2.0, 2.5));
        assert_eq!((comps3.n1, comps3.n2), (2.5, 2.5));
        assert_eq!((comps4.n1, comps4.n2), (2.5, 1.5));
        assert_eq!((comps5.n1, comps5.n2), (1.5, 1.0));
    }

    #[test]
    fn underpoint_is_offset_below_surface() {
        let ray = Ray::from_tuple((0., 0., -5.), (0., 0., 1.));
        let s = Shape::Sphere(Sphere::new(
            Some(Matrix4::<f32>::translation(0.0, 0.0, 1.0)),
            Some(Material::new_glass()),
            None,
            None,
        ));
        let i1 = Intersection::new(5.0, &s);
        let comps = i1.prepare_computation(&ray, None);
        assert_eq!(comps.under_point.z > EPSILON / 2.0, true);
        assert_eq!(comps.point.z < comps.under_point.z, true);
    }

    //P.161
    #[test]
    fn the_schlick_approximation_under_total_internal_reflection() {
        let shape = Shape::Sphere(Sphere::new(None, Some(Material::new_glass()), None, None));
        let ray = Ray::from_tuple((0., 0., 2_f32.sqrt() / 2.0), (0., 1., 0.));
        let i1 = Intersection::new(-2_f32.sqrt() / 2.0, &shape);
        let i2 = Intersection::new(2_f32.sqrt() / 2.0, &shape);
        let xs = vec![i1, i2];
        let comps = i2.prepare_computation(&ray, Some(&xs));
        let reflectance = comps.schlick();
        assert_relative_eq!(reflectance, 1.0);
    }

    #[test]
    fn the_schlick_approximation_with_a_perpendicular_viewing_angle() {
        let shape = Shape::Sphere(Sphere::new(None, Some(Material::new_glass()), None, None));
        let ray = Ray::from_tuple((0., 0., 0.), (0., 1., 0.));
        let i1 = Intersection::new(-1.0, &shape);
        let i2 = Intersection::new(1.0, &shape);
        let xs = vec![i1, i2];
        let comps = i2.prepare_computation(&ray, Some(&xs));
        let reflectance = comps.schlick();

        assert_relative_eq!(reflectance, 0.04);
    }

    #[test]
    fn shlick_approximation_with_small_angle() {
        let s = Shape::Sphere(Sphere::new(None, Some(Material::new_glass()), None, None));
        let ray = Ray::from_tuple((0., 0.99, -2.), (0., 0., 1.));
        let i1 = Intersection::new(1.8589, &s);
        let comps = i1.prepare_computation(&ray, None);
        let refrectance = comps.schlick();
        assert_relative_eq!(refrectance, 0.48873, epsilon = EPSILON);
    }
}
