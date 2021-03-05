extern crate approx;
extern crate nalgebra as na;
extern crate ulid;
use std::f32::consts::PI;
use std::iter::FromIterator;

use na::{Matrix4, Vector, Vector4};
use raytracer::{
    camera::Camera,
    canvas::Canvas,
    color::Color,
    cone::Cone,
    cube::Cube,
    cylinder::Cylinder,
    group::{Group, Scene},
    intersection::*,
    material::Material,
    matrix::CGMatrix,
    pattern::{Checker, Pattern, Ring, Stripe},
    plane::Plane,
    pointlight::PointLight,
    ppm::Ppm,
    ray::Ray,
    shape::Shape,
    sphere::Sphere,
    tuple::TupleOperation,
    world::World,
};
use vec_tree::VecTree;
fn main() {
    ch5();
    ch6();
    ch7();
    ch9();
    ch10();
    ch11();
    ch12();
    ch13_cylinder();
    ch13_cone();
    ch14();
}

fn ch5() {
    let ray_origin_z = -5.0;
    let wall_z = 10.0;
    let wall_size = 7.0;
    let canvas_pixels = 100.0;
    let pixel_size = wall_size / canvas_pixels;
    let half = wall_size / 2.0;

    let ray = Ray::from_tuple((0.0, 0.0, ray_origin_z), (0.0, 0.0, wall_z));
    let sphere = Sphere::new(Some(Matrix4::translation(0.1, 0.1, 0.0)), None, None, None);
    let shape = Shape::Sphere(sphere);

    let mut pixels: Vec<Color> = vec![];

    for y in 0..canvas_pixels as i32 {
        let world_y = half - pixel_size * (y as f32);
        for x in 0..canvas_pixels as i32 {
            let world_x = -half + pixel_size * (x as f32);
            let position = Vector4::point(world_x, world_y, wall_z);
            let direction = Vector::normalize(&(position - ray.origin));
            let r = Ray::new(ray.origin, direction);
            let intersections = shape.intersect(&r);
            let hit = intersections.hit();
            if hit.is_some() {
                let red = Color::rgb(1.0, 0.0, 0.0);
                pixels.push(red);
            } else {
                let black = Color::rgb(0.0, 0.0, 0.0);
                pixels.push(black);
            }
        }
    }
    let ppm = Ppm::new(canvas_pixels as usize, canvas_pixels as usize, &pixels);
    ppm.out("ch5");
}

fn ch6() {
    let material = Material::new(
        Some(Color::rgb(1.0, 0.2, 1.0)),
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    );
    let shape = Shape::Sphere(Sphere::new(
        Some(Matrix4::translation(0.1, 0.1, 0.0)),
        Some(material),
        None,
        None,
    ));
    let canvas = Canvas::render_single_shape(shape, 300);
    let ppm = Ppm::from_canvas(canvas);
    ppm.out("ch6");
}
fn ch7() {
    let floor_material = Material::new(
        Some(Color::rgb(1.0, 0.9, 0.9)),
        None,
        None,
        Some(0.0),
        None,
        None,
        None,
        None,
        None,
    );

    let floor = Shape::Sphere(Sphere::new(
        Some(Matrix4::scaling(10.0, 0.01, 10.0)),
        Some(floor_material.clone()),
        None,
        None,
    ));
    let left_wall = Shape::Sphere(Sphere::new(
        Some(
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(-PI / 4.0)
                * Matrix4::rotation_x(PI / 2.0)
                * Matrix4::scaling(10.0, 0.01, 10.0),
        ),
        Some(floor_material.clone()),
        None,
        None,
    ));
    let right_wall = Shape::Sphere(Sphere::new(
        Some(
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(PI / 4.0)
                * Matrix4::rotation_x(PI / 2.0)
                * Matrix4::scaling(10.0, 0.01, 10.0),
        ),
        Some(floor_material.clone()),
        None,
        None,
    ));

    let middle_sphere = Shape::Sphere(Sphere::new(
        Some(Matrix4::translation(-0.5, 1.0, 0.5)),
        Some(Material::new(
            Some(Color::rgb(1.0, 0.2, 0.2)),
            Some(0.7),
            Some(0.3),
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
    let right_sphere = Shape::Sphere(Sphere::new(
        Some(Matrix4::translation(1.5, 0.5, -0.5) * Matrix4::scaling(0.5, 0.5, 0.5)),
        Some(Material::new(
            Some(Color::rgb(0.3, 0.4, 1.0)),
            Some(0.7),
            Some(0.3),
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

    let world = World::new(
        Some(PointLight::new(
            Vector4::point(-10.0, 10.0, -10.0),
            Color::rgb(1.0, 1.0, 1.0),
        )),
        vec![floor, left_wall, right_wall, middle_sphere, right_sphere],
    );

    let camera = Camera::new(
        200.0,
        150.0,
        PI / 3.0,
        Vector4::point(0.0, 1.5, -5.0).view_transformation(
            &Vector4::point(0.0, 1.0, 0.0),
            &Vector4::vector(0.0, 1.0, 0.0),
        ),
    );

    let canvas = camera.render(&world);

    let ppm = canvas.to_ppm();
    ppm.out("ch7");
}

fn ch9() {
    let floor_material = Material::new(
        Some(Color::rgb(1.0, 0.9, 0.9)),
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    );

    let floor = Shape::Plane(Plane::new(None, Some(floor_material.clone()), None));

    let middle_sphere = Shape::Sphere(Sphere::new(
        Some(Matrix4::translation(-0.5, 1.0, 0.5)),
        Some(Material::new(
            Some(Color::rgb(1.0, 0.2, 0.2)),
            None,
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
    let right_sphere = Shape::Sphere(Sphere::new(
        Some(Matrix4::translation(1.5, 0.5, -0.5) * Matrix4::scaling(0.5, 0.5, 0.5)),
        Some(Material::new(
            Some(Color::rgb(0.3, 0.4, 1.0)),
            None,
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

    let world = World::new(
        Some(PointLight::new(
            Vector4::point(-10.0, 10.0, -10.0),
            Color::rgb(1.0, 1.0, 1.0),
        )),
        vec![floor, middle_sphere, right_sphere],
    );

    let ratio = 1.0;
    let camera = Camera::new(
        200.0 * ratio,
        150.0 * ratio,
        PI / 3.0,
        Vector4::point(0.0, 1.5, -5.0).view_transformation(
            &Vector4::point(0.0, 1.0, 0.0),
            &Vector4::vector(0.0, 1.0, 0.0),
        ),
    );

    let canvas = camera.render(&world);

    let ppm = canvas.to_ppm();
    ppm.out("ch9");
}

fn ch10() {
    let floor_material = Material::new(
        Some(Color::rgb(1.0, 0.9, 0.9)),
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        Some(Pattern::Checker(Checker::new(
            (1.0, 1.0, 1.0),
            (0.0, 0.0, 0.0),
            Matrix4::<f32>::scaling(0.5, 0.5, 0.5),
        ))),
    );

    let floor = Shape::Sphere(Sphere::new(
        Some(Matrix4::scaling(10.0, 0.01, 10.0)),
        Some(floor_material.clone()),
        None,
        None,
    ));

    let middle_sphere = Shape::Sphere(Sphere::new(
        Some(Matrix4::translation(-0.5, 1.0, 0.5)),
        Some(Material::new(
            Some(Color::rgb(1.0, 0.2, 0.2)),
            None,
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

    let right_back_sphere = Shape::Sphere(Sphere::new(
        Some(Matrix4::translation(1.7, 0.5, 2.5) * Matrix4::scaling(1.3, 1.3, 1.3)),
        Some(Material::new(
            Some(Color::rgb(0.3, 0.4, 1.0)),
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            Some(Pattern::Ring(Ring::new(
                (1.0, 0.0, 0.0),
                (1.0, 1.0, 1.0),
                Matrix4::<f32>::scaling(0.1, 0.1, 0.1),
            ))),
        )),
        None,
        None,
    ));

    let right_forward_sphere = Shape::Sphere(Sphere::new(
        Some(Matrix4::translation(1.2, 0.5, -0.9) * Matrix4::scaling(0.7, 0.7, 0.7)),
        Some(Material::new(
            Some(Color::rgb(0.0, 0.5, 1.0)),
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            Some(Pattern::Stripe(Stripe::new(
                (0.0, 0.8, 0.8),
                (1.0, 1.0, 1.0),
                Matrix4::<f32>::scaling(0.2, 0.2, 0.2),
            ))),
        )),
        None,
        None,
    ));

    let world = World::new(
        Some(PointLight::new(
            Vector4::point(-10.0, 10.0, -10.0),
            Color::rgb(1.0, 1.0, 1.0),
        )),
        vec![
            floor,
            middle_sphere,
            right_forward_sphere,
            right_back_sphere,
        ],
    );

    let ratio = 2.0;

    let camera = Camera::new(
        200.0 * ratio,
        150.0 * ratio,
        PI / 3.0,
        Vector4::point(0.0, 1.5, -5.0).view_transformation(
            &Vector4::point(0.0, 1.0, 0.0),
            &Vector4::vector(0.0, 1.0, 0.0),
        ),
    );

    let canvas = camera.render(&world);

    let ppm = canvas.to_ppm();
    ppm.out("ch10");
}

fn ch11() {
    let floor_material = Material::new(
        Some(Color::rgb(1.0, 0.9, 0.9)),
        None,
        None,
        Some(0.0),
        Some(0.0),
        Some(0.3),
        None,
        None,
        Some(Pattern::Checker(Checker::new(
            (1.0, 1.0, 1.0),
            (0.0, 0.0, 0.0),
            Matrix4::<f32>::scaling(0.5, 0.5, 0.5),
        ))),
    );

    let floor = Shape::Sphere(Sphere::new(
        Some(Matrix4::scaling(10.0, 0.01, 10.0)),
        Some(floor_material.clone()),
        None,
        None,
    ));

    let left_wall = Shape::Plane(Plane::new(
        Some(
            Matrix4::translation(0.0, 0.0, 8.0)
                * Matrix4::rotation_y(-PI / 4.0)
                * Matrix4::rotation_x(PI / 2.0),
        ),
        Some(Material::new(
            Some(Color::rgb(0.8, 0.8, 0.8)),
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
    let right_wall = Shape::Plane(Plane::new(
        Some(
            Matrix4::translation(0.0, 0.0, 8.0)
                * Matrix4::rotation_y(PI / 4.0)
                * Matrix4::rotation_x(PI / 2.0),
        ),
        Some(Material::new(
            Some(Color::rgb(0.8, 0.8, 0.8)),
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

    let middle_sphere = Shape::Sphere(Sphere::new(
        Some(Matrix4::translation(-0.5, 1.0, 0.5)),
        Some(Material::new(
            Some(Color::rgb(1.0, 0.2, 0.2)),
            None,
            None,
            None,
            None,
            Some(0.3),
            None,
            None,
            None,
        )),
        None,
        None,
    ));
    let right_back_sphere = Shape::Sphere(Sphere::new(
        Some(Matrix4::translation(1.5, 0.5, 1.5) * Matrix4::scaling(0.8, 0.8, 0.8)),
        Some(Material::new(
            Some(Color::rgb(0.3, 0.4, 1.0)),
            None,
            None,
            None,
            None,
            None,
            None,
            Some(0.8),
            Some(Pattern::Stripe(Stripe::new(
                (0.0, 0.3, 1.0),
                (1.0, 1.0, 1.0),
                Matrix4::<f32>::scaling(0.7, 0.7, 0.7),
            ))),
        )),
        None,
        None,
    ));
    let right_forward_sphere = Shape::Sphere(Sphere::new(
        Some(Matrix4::translation(1.5, 0.5, -0.5) * Matrix4::scaling(0.3, 0.3, 0.3)),
        Some(Material::new(
            Some(Color::rgb(0.1, 0.9, 0.2)),
            None,
            None,
            None,
            None,
            Some(0.4),
            None,
            None,
            None,
        )),
        None,
        None,
    ));

    let world = World::new(
        Some(PointLight::new(
            Vector4::point(-10.0, 10.0, -10.0),
            Color::rgb(1.0, 1.0, 1.0),
        )),
        vec![
            floor,
            left_wall,
            right_wall,
            middle_sphere,
            right_forward_sphere,
            right_back_sphere,
        ],
    );

    let ratio = 3.0;
    let camera = Camera::new(
        200.0 * ratio,
        150.0 * ratio,
        PI / 3.0,
        Vector4::point(0.0, 1.5, -5.0).view_transformation(
            &Vector4::point(0.0, 1.0, 0.0),
            &Vector4::vector(0.0, 1.0, 0.0),
        ),
    );

    let canvas = camera.render(&world);

    let ppm = canvas.to_ppm();
    ppm.out("ch11");
}

fn ch12() {
    let floor_material = Material::new(
        Some(Color::rgb(1.0, 0.9, 0.9)),
        None,
        None,
        Some(0.0),
        Some(0.0),
        Some(0.7),
        Some(0.0),
        Some(0.0),
        Some(Pattern::Checker(Checker::new(
            (1.0, 1.0, 1.0),
            (0.0, 0.0, 0.0),
            Matrix4::<f32>::scaling(0.5, 0.5, 0.5),
        ))),
    );

    let floor = Shape::Sphere(Sphere::new(
        Some(Matrix4::scaling(10.0, 0.01, 10.0)),
        Some(floor_material.clone()),
        None,
        None,
    ));

    let shapes_iter = [-7., -5., -3., -1., 1., 3., 5., 7.].iter().map(|x| {
        if *x < 0. {
            Shape::Cube(Cube::new(
                Some(
                    Matrix4::translation(0.45 * x, 0.3, 0.6 * x.abs() - 3.0)
                        * Matrix4::scaling(0.3, 0.3, 0.3)
                        * Matrix4::rotation_y(0.45),
                ),
                Some(Material::new(
                    Some(Color::rgb(
                        1.5 * x.abs() / 12.,
                        0.1,
                        1.0 - x.abs() * 1.5 / 12.,
                    )),
                    None,
                    None,
                    None,
                    None,
                    Some(0.4),
                    None,
                    None,
                    None,
                )),
            ))
        } else {
            Shape::Sphere(Sphere::new(
                Some(
                    Matrix4::translation(0.45 * x, 0.4, 0.6 * x.abs() - 3.0)
                        * Matrix4::scaling(0.4, 0.4, 0.4)
                        * Matrix4::rotation_y(0.45),
                ),
                Some(Material::new(
                    Some(Color::rgb(
                        1.5 * x.abs() / 12.,
                        0.1,
                        1.0 - x.abs() * 1.5 / 12.,
                    )),
                    None,
                    None,
                    None,
                    None,
                    Some(0.4),
                    None,
                    None,
                    None,
                )),
                None,
                None,
            ))
        }
    });
    let mut shapes = Vec::from_iter(shapes_iter);
    shapes.push(floor);
    let world = World::new(
        Some(PointLight::new(
            Vector4::point(-10.0, 10.0, -10.0),
            Color::rgb(1.0, 1.0, 1.0),
        )),
        shapes,
    );

    let ratio = 1.0;
    let camera = Camera::new(
        200.0 * ratio,
        150.0 * ratio,
        PI / 3.0,
        Vector4::point(0.0, 1.0, -5.0).view_transformation(
            &Vector4::point(0.0, 0.5, 0.0),
            &Vector4::vector(0.0, 0.5, 0.0),
        ),
    );

    let canvas = camera.render(&world);

    let ppm = canvas.to_ppm();
    ppm.out("ch12");
}

fn ch13_cylinder() {
    let floor_material = Material::new(
        Some(Color::rgb(1.0, 0.9, 0.9)),
        None,
        None,
        Some(0.0),
        Some(0.0),
        Some(0.3),
        Some(0.0),
        Some(0.0),
        Some(Pattern::Checker(Checker::new(
            (1.0, 1.0, 1.0),
            (0.0, 0.0, 0.0),
            Matrix4::<f32>::scaling(0.5, 0.5, 0.5),
        ))),
    );

    let floor = Shape::Sphere(Sphere::new(
        Some(Matrix4::scaling(10.0, 0.01, 10.0)),
        Some(floor_material.clone()),
        None,
        None,
    ));
    let shapes_iter = [-5., -3., -1., 1., 3., 5.].iter().map(|x| {
        Shape::Cylinder(Cylinder::new(
            Some(
                Matrix4::translation(0.35 * x, 0.5, -0.4 * x.abs())
                    * Matrix4::scaling(0.3, 1.0, 0.3),
            ),
            Some(Material::new(
                Some(Color::rgb(0.1, x.abs() / 12., 0.9)),
                None,
                None,
                None,
                None,
                Some(0.4),
                None,
                None,
                None,
            )),
            Some(0.0),
            Some(1.5),
            Some(true),
        ))
    });
    let mut shapes = Vec::from_iter(shapes_iter);
    shapes.push(floor);
    let world = World::new(
        Some(PointLight::new(
            Vector4::point(-10.0, 10.0, -10.0),
            Color::rgb(1.0, 1.0, 1.0),
        )),
        shapes,
    );

    let ratio = 3.0;
    let camera = Camera::new(
        200.0 * ratio,
        150.0 * ratio,
        PI / 3.0,
        Vector4::point(0.0, 1.5, -5.0).view_transformation(
            &Vector4::point(0.0, 1.0, 0.0),
            &Vector4::vector(0.0, 1.0, 0.0),
        ),
    );

    let canvas = camera.render(&world);

    let ppm = canvas.to_ppm();
    ppm.out("ch13_cylinder");
}

fn ch13_cone() {
    let floor_material = Material::new(
        Some(Color::rgb(1.0, 0.9, 0.9)),
        None,
        None,
        Some(0.0),
        Some(0.0),
        Some(0.3),
        Some(0.0),
        Some(0.0),
        Some(Pattern::Checker(Checker::new(
            (1.0, 1.0, 1.0),
            (0.0, 0.0, 0.0),
            Matrix4::<f32>::scaling(0.5, 0.5, 0.5),
        ))),
    );

    let floor = Shape::Sphere(Sphere::new(
        Some(Matrix4::scaling(10.0, 0.01, 10.0)),
        Some(floor_material.clone()),
        None,
        None,
    ));
    let shapes_iter = [-5., -3., -1., 1., 3., 5.].iter().map(|x| {
        Shape::Cone(Cone::new(
            Some(
                Matrix4::translation(0.35 * x, 1.0, -0.4 * x.abs())
                    * Matrix4::scaling(0.4, 1.0, 0.4),
            ),
            Some(Material::new(
                Some(Color::rgb(0.1, x.abs() / 12., 0.9)),
                None,
                None,
                None,
                None,
                Some(0.4),
                None,
                None,
                None,
            )),
            Some(0.0),
            Some(1.0),
            Some(true),
        ))
    });
    let mut shapes = Vec::from_iter(shapes_iter);
    shapes.push(floor);
    let world = World::new(
        Some(PointLight::new(
            Vector4::point(-10.0, 10.0, -10.0),
            Color::rgb(1.0, 1.0, 1.0),
        )),
        shapes,
    );

    let ratio = 3.0;
    let camera = Camera::new(
        200.0 * ratio,
        150.0 * ratio,
        PI / 3.0,
        Vector4::point(0.0, 1.0, -5.0).view_transformation(
            &Vector4::point(0.0, 1.0, 0.0),
            &Vector4::vector(0.0, 1.0, 0.0),
        ),
    );

    let canvas = camera.render(&world);

    let ppm = canvas.to_ppm();
    ppm.out("ch13_cone");
}

fn ch14() {
    let floor_material = Material::new(
        Some(Color::rgb(1.0, 0.9, 0.9)),
        None,
        None,
        Some(0.0),
        Some(0.0),
        Some(0.7),
        Some(0.0),
        Some(0.0),
        Some(Pattern::Checker(Checker::new(
            (1.0, 1.0, 1.0),
            (0.0, 0.0, 0.0),
            Matrix4::<f32>::scaling(0.5, 0.5, 0.5),
        ))),
    );

    let floor = Shape::Sphere(Sphere::new(
        Some(Matrix4::scaling(10.0, 0.01, 10.0)),
        Some(floor_material.clone()),
        None,
        None,
    ));

    let c1 = Shape::Cube(Cube::new(
        Some(Matrix4::translation(-3.0, 1.5, 1.0)),
        Some(Material::new(
            None,
            None,
            None,
            None,
            None,
            Some(0.4),
            None,
            None,
            None,
        )),
    ));

    let c2 = Shape::Cube(Cube::new(
        Some(Matrix4::translation(1.0, 1.5, 1.0)),
        Some(Material::new(
            Some(Color::rgb(0.8, 0.1, 0.1)),
            None,
            None,
            None,
            None,
            Some(0.4),
            None,
            None,
            None,
        )),
    ));
    let c3 = Shape::Cube(Cube::new(
        Some(Matrix4::translation(2.0, 1., 0.5)),
        Some(Material::new(
            Some(Color::rgb(0.1, 0.1, 0.9)),
            None,
            None,
            None,
            None,
            Some(0.4),
            None,
            None,
            None,
        )),
    ));

    let g = Shape::Group(Group::new(
        Some(
            Matrix4::translation(0.3, 1.0, 1.0)
                * Matrix4::scaling(0.4, 0.4, 0.4)
                * Matrix4::rotation_y(0.45)
                * Matrix4::rotation_x(0.45)
                * Matrix4::rotation_z(0.45),
        ),
        None,
    ));

    let mut tree = VecTree::new();
    let fln = tree.insert_root(floor.id());
    let gn = tree.insert(g.id(), fln);
    let c1n = tree.insert(c1.id(), gn);
    let c2n = tree.insert(c2.id(), gn);
    let c3n = tree.insert(c3.id(), fln);

    let scene = Scene::from_shapes(
        tree,
        vec![floor.clone(), g.clone(), c1.clone(), c2.clone(), c3.clone()],
    );
    let shapes = scene.to_transformed_shapes();

    let world = World::new(
        Some(PointLight::new(
            Vector4::point(-10.0, 10.0, -10.0),
            Color::rgb(1.0, 1.0, 1.0),
        )),
        shapes,
    );

    let ratio = 1.0;
    let camera = Camera::new(
        200.0 * ratio,
        150.0 * ratio,
        PI / 3.0,
        Vector4::point(0.0, 1.0, -5.0).view_transformation(
            &Vector4::point(0.0, 0.5, 0.0),
            &Vector4::vector(0.0, 0.5, 0.0),
        ),
    );

    let canvas = camera.render(&world);

    let ppm = canvas.to_ppm();
    ppm.out("ch14");
}
