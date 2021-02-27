use core::f32;

use computation::Computations;
use na::Vector4;

use crate::{
    color::Color,
    computation,
    epsilon::EPSILON,
    intersection::{Intersection, Intersections},
    pointlight::PointLight,
    ray::Ray,
    shape::Shape,
};

#[derive(Clone, Debug)]
pub struct World {
    pub pointlight: PointLight,
    pub shapes: Vec<Shape>,
}

impl World {
    pub fn new(pointlight: Option<PointLight>, shapes: Vec<Shape>) -> World {
        World {
            pointlight: match pointlight {
                Some(x) => x,
                None => PointLight::new_default(),
            },
            shapes,
        }
    }
    pub fn new_default() -> World {
        World::new(None, Shape::new_world_defaults())
    }
    pub fn intersect(&self, ray: &Ray) -> Vec<Intersection> {
        let mut intersections = self
            .shapes
            .iter()
            .map(|x| x.intersect(&ray))
            .flatten()
            .collect::<Vec<_>>();
        intersections.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
        intersections.clone()
    }
    pub fn is_shadowed(&self, point: Vector4<f32>) -> bool {
        let light_vector = self.pointlight.position - point;
        let distance = light_vector.magnitude();
        let direction = light_vector.normalize();
        let ray = Ray::new(point, direction);
        let intersections = self.intersect(&ray);
        let hit = intersections.hit();
        if hit.is_some() && hit.unwrap().t < distance + EPSILON {
            return true;
        } else {
            return false;
        }
    }
    pub fn color_at(&self, ray: &Ray, remaining: usize) -> Color {
        let intersections = self.intersect(&ray);
        let hit = intersections.hit();
        if hit.is_none() {
            return Color::rgb(0.0, 0.0, 0.0);
        }
        let computations = hit.unwrap().prepare_computation(&ray, Some(&intersections));
        self.shade_hit(&computations, remaining)
    }

    pub fn reflected_color(&self, computations: &Computations, remaining: usize) -> Color {
        if remaining <= 0 || computations.shape.material().reflective == 0.0 {
            return Color::rgb(0.0, 0.0, 0.0);
        }
        let reflected_ray = Ray::new(computations.over_point, computations.reflectv);
        let color = self.color_at(&reflected_ray, remaining - 1);
        let reflected_color = color * computations.shape.material().reflective;
        reflected_color
    }

    pub fn refracted_color(&self, computations: &Computations, remaining: usize) -> Color {
        if remaining <= 0 || computations.shape.material().transparency == 0.0 {
            return Color::rgb(0.0, 0.0, 0.0);
        }
        // Snell's law: sin(theta_i) / sin(theta_t) = eta_2 / eta_1gg
        // s.t. n2 := eta_2, n1 := eta_1
        // Following ratio is inversion of Snell's law
        let n_ratio = computations.n1 / computations.n2;
        // cos(theta_i)
        let cos_ti = computations.eyev.dot(&computations.normalv);

        // sin(theta_t) **
        let sin2_t = n_ratio.powi(2) * (1.0 - cos_ti.powi(2));
        if sin2_t > 1.0 {
            return Color::rgb(0.0, 0.0, 0.0);
        }
        let cos_t = (1.0 - sin2_t).sqrt();
        let direction =
            (n_ratio * cos_ti - cos_t) * computations.normalv - n_ratio * computations.eyev;

        let refrect_ray = Ray::new(computations.under_point, direction);
        let color =
            self.color_at(&refrect_ray, remaining - 1) * computations.shape.material().transparency;
        color
    }

    pub fn shade_hit(&self, computations: &Computations, remaining: usize) -> Color {
        let shadowed = self.is_shadowed(computations.over_point);
        let material = computations.shape.material();

        let surface_color = material.lighting(
            &self.pointlight,
            computations.over_point,
            computations.eyev,
            computations.normalv,
            shadowed,
            Some(&computations.shape),
        );

        let reflected_color = self.reflected_color(&computations, remaining);
        let refracted_color = self.refracted_color(&computations, remaining);
        let material = computations.shape.material();
        if material.reflective > 0.0 && material.transparency > 0.0 {
            let reflectance = computations.schlick();
            return surface_color
                + reflected_color * reflectance
                + refracted_color * (1. - reflectance);
        }
        surface_color + refracted_color + reflected_color
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tuple::TupleOperation;
    use crate::{
        epsilon,
        material::Material,
        matrix::CGMatrix,
        pattern::{Pattern, Test},
        plane::Plane,
        shape::Shape,
        sphere::Sphere,
    };
    use approx::assert_relative_eq;
    use na::Matrix4;

    #[test]
    fn intersect_world_with_ray() {
        let world = World::new_default();
        let ray = Ray::from_tuple((0., 0., -5.), (0., 0., 1.));
        let xs = world.intersect(&ray);
        assert_eq!(xs.len(), 4);
        assert_eq!(xs[0].t, 4.0);
        assert_eq!(xs[1].t, 4.5);
        assert_eq!(xs[2].t, 5.5);
        assert_eq!(xs[3].t, 6.0);
    }

    // P.95
    #[test]
    fn shadeing_intersection() {
        let world = World::new_default();
        let i = Intersection::new(4.0, &world.shapes[0]);
        let ray = Ray::from_tuple((0., 0., -5.), (0., 0., 1.));
        let computations = i.prepare_computation(&ray, None);
        let c = world.shade_hit(&computations, 2);
        assert_eq!(c, Color::rgb(0.38066, 0.47583, 0.2855));
    }

    // P.95
    #[test]
    fn shadeing_intersection_from_inside() {
        let mut world = World::new_default();
        world.pointlight = PointLight::from_tuple((0., 0.25, 0.), (1., 1., 1.));
        let i = Intersection::new(0.5, &world.shapes[1]);
        let ray = Ray::from_tuple((0., 0., 0.), (0., 0., 1.));
        let computations = i.prepare_computation(&ray, None);
        let c = world.shade_hit(&computations, 1);
        assert_eq!(c, Color::rgb(0.90498, 0.90498, 0.90498));
    }

    // P.114
    #[test]
    fn shadehit_is_given_intersection_in_shadow() {
        let s1 = Shape::Sphere(Sphere::new_default());
        let s2 = Shape::Sphere(Sphere::new(
            Some(Matrix4::<f32>::translation(0.0, 0.0, 10.0)),
            None,
            None,
            None,
            None,
        ));
        let world = World::new(
            Some(PointLight::from_tuple((0.0, 0.0, -10.0), (1.0, 1.0, 1.0))),
            vec![s1.clone(), s2.clone()],
        );
        let ray = Ray::from_tuple((0., 0., 5.), (0., 0., 1.));
        let i = Intersection::new(4.0, &s2);
        let computations = i.prepare_computation(&ray, None);
        let c = world.shade_hit(&computations, 2);
        assert_eq!(c, Color::rgb(0.1, 0.1, 0.1));
    }

    // P.115
    #[test]
    fn hit_offset_point() {
        let shape = Shape::Sphere(Sphere::new(
            Some(Matrix4::<f32>::translation(0.0, 0.0, 1.0)),
            None,
            None,
            None,
            None,
        ));
        let i = Intersection::new(5.0, &shape);
        let ray = Ray::from_tuple((0., 0., -5.), (0., 0., 1.));

        let computations = i.prepare_computation(&ray, None);
        assert_eq!(computations.over_point[2] < -epsilon::EPSILON / 2.0, true);
        assert_eq!(computations.point[2] > computations.over_point.z, true);
    }

    //P.155
    #[test]
    fn refracted_color_with_opaque_surface() {
        let world = World::new_default();
        let ray = Ray::from_tuple((0., 0., -5.0), (0., 0., 1.));
        let i1 = Intersection::new(4.0, &world.shapes[0]);
        let i2 = Intersection::new(6.0, &world.shapes[0]);
        let xs = vec![i1, i2];
        let comps = xs[0].prepare_computation(&ray, Some(&xs));
        let color = world.refracted_color(&comps, 5);
        assert_relative_eq!(color.vector(), Color::rgb(0., 0., 0.).vector());
    }

    //P.156
    #[test]
    fn refracted_color_under_total_internal_reflection() {
        let mut shapes = Shape::new_world_defaults();
        let mut a = shapes[0].clone();
        let mut am = a.material().clone();
        am.transparency = 1.0;
        am.refractive = 1.5;
        a.set_material(am);
        shapes.push(a.clone());

        let ray = Ray::from_tuple((0., 0., 2_f32.sqrt() / 2.0), (0., 1., 0.));
        let i1 = Intersection::new(-2_f32.sqrt() / 2.0, &a);
        let i2 = Intersection::new(2_f32.sqrt() / 2.0, &a);
        let xs = vec![i1, i2];
        let comps = xs[1].prepare_computation(&ray, Some(&xs));
        let world = World::new(None, shapes);
        let color = world.refracted_color(&comps, 5);
        assert_relative_eq!(
            color.vector(),
            Color::rgb(0., 0., 0.).vector(),
            epsilon = EPSILON
        );
    }

    //P.158
    #[test]
    fn refracted_color_refracted_ray() {
        let shapes = Shape::new_world_defaults();
        let mut a = shapes[0].clone();
        let mut am = a.material().clone();
        am.ambient = 1.0;
        am.pattern = Some(Pattern::Test(Test::new()));
        a.set_material(am);

        let mut b = shapes[1].clone();
        let mut bm = b.material().clone();
        bm.transparency = 1.0;
        bm.refractive = 1.5;
        b.set_material(bm);

        let ray = Ray::from_tuple((0., 0., 0.1), (0., 1., 0.));
        let i1 = Intersection::new(-0.9899, &a);
        let i2 = Intersection::new(-0.4899, &b);
        let i3 = Intersection::new(0.4899, &b);
        let i4 = Intersection::new(0.9899, &a);
        let xs = vec![i1, i2, i3, i4];
        let comps = xs[2].prepare_computation(&ray, Some(&xs));
        let world = World::new(None, vec![a.clone(), b.clone()]);
        let color = world.refracted_color(&comps, 5);
        assert_relative_eq!(
            color.vector(),
            Color::rgb(0., 0.99888, 0.04725).vector(),
            epsilon = EPSILON * 10., //un accurate?
        );
    }

    //P.159
    #[test]
    fn shade_hit_transparent_material() {
        let floor = Shape::Plane(Plane::new(
            Some(Matrix4::<f32>::translation(0., -1., 0.)),
            Some(Material::new(
                None,
                None,
                None,
                None,
                None,
                None,
                Some(1.5),
                Some(0.5),
                None,
            )),
            None,
            None,
        ));
        let ball = Shape::Sphere(Sphere::new(
            Some(Matrix4::<f32>::translation(0., -3.5, -0.5)),
            Some(Material::new(
                Some(Color::rgb(1.0, 0.0, 0.0)),
                Some(0.5),
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
            None,
        ));
        let mut shapes = Shape::new_world_defaults();
        shapes.push(floor.clone());
        shapes.push(ball.clone());
        assert_eq!(shapes.len(), 4);
        let world = World::new(None, shapes);
        let i1 = Intersection::new(2_f32.sqrt(), &floor);
        let xs = vec![i1];
        let ray = Ray::from_tuple(
            (0., 0., -3.),
            (0., -(2_f32.sqrt() / 2.0), 2_f32.sqrt() / 2.0),
        );
        let computations = xs[0].prepare_computation(&ray, Some(&xs));
        assert_relative_eq!(computations.n1, 1.0);
        assert_relative_eq!(
            computations.point,
            Vector4::point(0., -0.9999990000000002, -1.9999999999999998),
            epsilon = EPSILON
        );
        assert_relative_eq!(
            computations.eyev,
            Vector4::vector(0., 0.7071067811865476, -0.7071067811865476),
            epsilon = EPSILON
        );
        assert_relative_eq!(computations.normalv, Vector4::vector(0., 1., 0.));
        assert_relative_eq!(
            computations.reflectv,
            Vector4::vector(0., 0.7071067811865476, 0.7071067811865476),
            epsilon = EPSILON
        );
        assert_eq!(computations.inside, false);
        assert_relative_eq!(
            computations.under_point,
            Vector4::point(0., -1.0000010000000001, -1.9999999999999998),
            epsilon = EPSILON
        );
        assert_eq!(computations.n1, 1.0);

        let color = world.shade_hit(&computations, 5);

        assert_relative_eq!(
            color.vector(),
            Color::rgb(0.93642, 0.68642, 0.68642).vector(),
            epsilon = EPSILON
        );
    }

    //P.164
    #[test]
    fn shade_hit_reflective_transparent_material() {
        let ray = Ray::from_tuple(
            (0., 0., -3.),
            (0., -(2_f32.sqrt() / 2.0), 2_f32.sqrt() / 2.0),
        );
        let floor = Shape::Plane(Plane::new(
            Some(Matrix4::<f32>::translation(0., -1., 0.)),
            Some(Material::new(
                None,
                None,
                None,
                None,
                None,
                Some(0.5),
                Some(1.5),
                Some(0.5),
                None,
            )),
            None,
            None,
        ));
        let ball = Shape::Sphere(Sphere::new(
            Some(Matrix4::<f32>::translation(0., -3.5, -0.5)),
            Some(Material::new(
                Some(Color::rgb(1.0, 0.0, 0.0)),
                Some(0.5),
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
            None,
        ));
        let world = World::new(None, vec![floor.clone(), ball.clone()]);
        let xs = Intersection::new(2_f32.sqrt(), &floor);
        let computations = xs.prepare_computation(&ray, Some(&vec![xs]));
        let color = world.shade_hit(&computations, 5);
        assert_relative_eq!(
            color.vector(),
            Color::rgb(0.93391, 0.69643, 0.69243).vector(),
            epsilon = EPSILON * 100.0
        );
    }
}
