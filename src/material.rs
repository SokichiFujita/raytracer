use crate::tuple::TupleOperation;
use crate::{color::Color, refractive::RefractiveIndex, shape::Shape};
use crate::{pattern::Pattern, pointlight::PointLight};
use na::Vector4;

#[derive(Copy, Clone, PartialEq, Debug)]
pub struct Material {
    pub color: Color,
    pub ambient: f32,
    pub diffuse: f32,
    pub specular: f32,
    pub shininess: f32,
    pub reflective: f32,
    pub refractive: f32,
    pub transparency: f32,
    pub pattern: Option<Pattern>,
}

impl Material {
    pub fn new(
        color: Option<Color>,
        ambient: Option<f32>,
        diffuse: Option<f32>,
        specular: Option<f32>,
        shininess: Option<f32>,
        reflective: Option<f32>,
        refractive: Option<f32>,
        transparency: Option<f32>,
        pattern: Option<Pattern>,
    ) -> Material {
        Material {
            color: if color == None {
                Color::rgb(1.0, 1.0, 1.0)
            } else {
                color.unwrap()
            },
            ambient: if ambient == None {
                0.1
            } else {
                ambient.unwrap()
            },
            diffuse: if diffuse == None {
                0.9
            } else {
                diffuse.unwrap()
            },
            specular: if specular == None {
                0.9
            } else {
                specular.unwrap()
            },
            shininess: if shininess == None {
                200.0
            } else {
                shininess.unwrap()
            },
            reflective: if reflective == None {
                0.0
            } else {
                reflective.unwrap()
            },
            refractive: if refractive == None {
                RefractiveIndex::VACCUME
            } else {
                refractive.unwrap()
            },
            transparency: if transparency == None {
                0.0
            } else {
                transparency.unwrap()
            },
            pattern: pattern,
        }
    }

    pub fn new_default() -> Material {
        Material::new(None, None, None, None, None, None, None, None, None)
    }

    pub fn new_default_world() -> Material {
        Material::new(
            Some(Color::rgb(0.8, 1.0, 0.6)),
            None,
            Some(0.7),
            Some(0.2),
            None,
            None,
            None,
            None,
            None,
        )
    }

    pub fn new_glass() -> Material {
        Material::new(
            None,
            None,
            None,
            None,
            None,
            None,
            Some(1.5),
            Some(1.0),
            None,
        )
    }

    pub fn lighting(
        &self,
        light: &PointLight,
        point: Vector4<f32>,
        eyev: Vector4<f32>,
        normalv: Vector4<f32>,
        is_shadow: bool,
        shape: Option<&Shape>,
    ) -> Color {
        let c = match &self.pattern {
            Some(x) => x.pattern_at(shape.unwrap(), point),
            None => self.color.clone(),
        };
        // Add light intensity to material surface
        let effective_color = c.hadamard_product(&light.intensity);

        // Calcurate ambient lighting which is reflected light or environment light
        let ambient = (effective_color * self.ambient).vector();
        if is_shadow {
            return Color::rgb(ambient.x, ambient.y, ambient.z);
        }

        // Direction from light origin to given point
        let lightv = (light.position - point).normalize();

        // Dot product means cos of given two vectors. It represents angle of two vectors
        // ex.) two unit vectors are
        //  0 degree (same) => 1
        //  45 degree => √2/2
        //  90 degree => 0
        //  180 degree (opposite directions) => -1
        let light_dot_normal = lightv.dot(&normalv);

        let mut diffuse = Color::rgb(0.0, 0.0, 0.0).vector();
        let mut specular = Color::rgb(0.0, 0.0, 0.0).vector();

        // Angle of light and normal less than or equal to 90 degree
        if light_dot_normal >= 0.0 {
            // if light and normal is close then color will not be dark
            // if light and normal is far then color will be dark
            let k = self.diffuse * light_dot_normal;
            diffuse = k * effective_color.vector();

            // Calculate how reflecting vector and eye vector are close
            let reflectv = (-1.0 * lightv).reflect(normalv);
            let reflect_dot_eye = reflectv.dot(&eyev);

            // Angle of reflecting vector and eye vector are less than 90 degree
            if reflect_dot_eye > 0.0 {
                let factor = reflect_dot_eye.powf(self.shininess);
                specular = factor * self.specular * light.intensity.vector();
            }
        }
        Color::from_vector(ambient + diffuse + specular)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{sphere::Sphere, world::World};

    #[test]
    fn default_material() {
        let m = Material::new_default();
        assert_eq!(m.color, Color::rgb(1., 1., 1.));
        assert_eq!(m.ambient, 0.1);
        assert_eq!(m.diffuse, 0.9);
        assert_eq!(m.specular, 0.9);
        assert_eq!(m.shininess, 200.0);
    }

    #[test]
    fn sphere_has_default_material() {
        let s = Sphere::new_default();
        let m = Material::new_default();
        assert_eq!(s.material, m);
    }

    #[test]
    fn sphere_may_be_assigned_a_material() {
        let mut s = Sphere::new_default();
        let mut m = Material::new_default();
        m.ambient = 1.0;
        s.material = m.clone();
        assert_eq!(s.material, m);
    }

    //Ch6 P86
    #[test]
    fn lighting_light_0_eye_0() {
        let m = Material::new_default();
        let position = Vector4::point(0., 0., 0.);
        let eyev = Vector4::vector(0., 0., -1.);
        let normalv = Vector4::vector(0., 0., -1.);
        let light = PointLight::from_tuple((0., 0., -10.), (1., 1., 1.));
        let result = m.lighting(&light, position, eyev, normalv, false, None);
        assert_eq!(result, Color::rgb(1.9, 1.9, 1.9))
    }

    //Ch6 P86
    #[test]
    fn lighting_light_0_eye_45() {
        let m = Material::new_default();
        let position = Vector4::point(0., 0., 0.);
        let eyev = Vector4::vector(0., 2_f32.sqrt() / 2., 2_f32.sqrt() / 2.);
        let normalv = Vector4::vector(0., 0., -1.);
        let light = PointLight::from_tuple((0., 0., -10.), (1., 1., 1.));
        let result = m.lighting(&light, position, eyev, normalv, false, None);
        assert_eq!(result, Color::rgb(1., 1., 1.))
    }

    //Ch6 P87
    #[test]
    fn lighting_light_45_eye_0() {
        let m = Material::new_default();
        let position = Vector4::point(0., 0., 0.);
        let eyev = Vector4::vector(0., 0., -1.);
        let normalv = Vector4::vector(0., 0., -1.);
        let light = PointLight::from_tuple((0., 10., -10.), (1., 1., 1.));
        let result = m.lighting(&light, position, eyev, normalv, false, None);
        assert_eq!(result, Color::rgb(0.7364, 0.7364, 0.7364));
    }

    //Ch6 P87
    #[test]
    fn lighting_light_45_eye_315() {
        let m = Material::new_default();
        let position = Vector4::point(0., 0., 0.);
        let eyev = Vector4::vector(0., -2_f32.sqrt() / 2., -2_f32.sqrt() / 2.);
        let normalv = Vector4::vector(0., 0., -1.);
        let light = PointLight::from_tuple((0., 10., -10.), (1., 1., 1.));
        let result = m.lighting(&light, position, eyev, normalv, false, None);
        assert_eq!(result, Color::rgb(1.6364, 1.6364, 1.6364));
    }

    //Ch6 P88
    #[test]
    fn lighting_light_180_eye_0() {
        let m = Material::new_default();
        let position = Vector4::point(0., 0., 0.);
        let eyev = Vector4::vector(0., 0., -1.);
        let normalv = Vector4::vector(0., 0., -1.);
        let light = PointLight::from_tuple((0., 0., 10.), (1., 1., 1.));
        let result = m.lighting(&light, position, eyev, normalv, false, None);
        assert_eq!(result, Color::rgb(0.1, 0.1, 0.1));
    }

    //Ch8 P110
    #[test]
    fn lighting_with_surface_in_shadow() {
        let m = Material::new_default();
        let position = Vector4::point(0., 0., 0.);
        let eyev = Vector4::vector(0., 0., -1.);
        let normalv = Vector4::vector(0., 0., -1.);
        let light = PointLight::from_tuple((0., 0., -10.), (1., 1., 1.));
        let result = m.lighting(&light, position, eyev, normalv, true, None);
        assert_eq!(result, Color::rgb(0.1, 0.1, 0.1));
    }

    //Ch8 P111
    #[test]
    fn no_shadow_nothing_collinear_point_light() {
        let w = World::new_default();
        let p = Vector4::point(0., 10., 0.);
        assert_eq!(w.is_shadowed(p), false);
    }

    //Ch8 P111
    #[test]
    fn shadow_object_between_point_and_light() {
        let w = World::new_default();
        let p = Vector4::point(10., -10., 10.);
        assert_eq!(w.is_shadowed(p), true);
    }

    //Ch8 P111
    #[test]
    fn no_shadow_object_behind_light() {
        let w = World::new_default();
        let p = Vector4::point(-20., 20., -20.);
        assert_eq!(w.is_shadowed(p), false);
    }

    //Ch8 P111
    #[test]
    fn no_shadow_object_behind_point() {
        let w = World::new_default();
        let p = Vector4::point(-2., 2., -2.);
        assert_eq!(w.is_shadowed(p), false);
    }

    //P150
    #[test]
    fn transparency_refractive_index_default_material() {
        let m = Material::new_default();
        assert_eq!(m.transparency, 0.0);
        assert_eq!(m.refractive, 1.0);
    }

    //P151
    #[test]
    fn helper_glass_material() {
        let m = Material::new_glass();
        assert_eq!(m.transparency, 1.0);
        assert_eq!(m.refractive, 1.5);
    }
}
