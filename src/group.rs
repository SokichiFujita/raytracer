use std::collections::HashMap;

use crate::{intersection::Intersection, material::Material, ray::Ray, shape::Shape};
use na::Matrix4;
use nalgebra::Vector4;
use vec_tree::{Index, VecTree};

#[derive(Clone, Debug)]
pub struct Scene {
    tree: VecTree<String>,
    shapes: HashMap<String, Shape>,
}

impl Scene {
    pub fn new(tree: VecTree<String>, shapes: HashMap<String, Shape>) -> Scene {
        Scene {
            tree: tree,
            shapes: shapes,
        }
    }

    pub fn new_default() -> Scene {
        Scene::new(VecTree::new(), HashMap::new())
    }

    pub fn get_shape_by_id(&self, id: &String) -> Option<&Shape> {
        self.shapes.get(id)
    }

    pub fn get_shape_by_index(&self, index: Index) -> Option<&Shape> {
        match self.tree.get(index) {
            Some(x) => self.shapes.get(x),
            None => None,
        }
    }

    pub fn from_shapes(tree: VecTree<String>, shapes: Vec<Shape>) -> Scene {
        let mut shapes_hashmap = HashMap::new();
        shapes.iter().for_each(|x| {
            shapes_hashmap.insert(x.id(), x.clone());
        });
        let scene = Scene::new(tree.clone(), shapes_hashmap);
        scene
    }

    pub fn len_children(&self, index: Index) -> usize {
        let len = self
            .tree
            .children(index)
            .map(|node| self.tree.get(node).unwrap())
            .collect::<Vec<&String>>()
            .len();
        len
    }
    pub fn children(&self, index: Index) -> Vec<&String> {
        let children = self
            .tree
            .children(index)
            .map(|node| self.tree.get(node).unwrap())
            .collect::<Vec<&String>>();
        children
    }
    pub fn parent(&self, index: Index) -> Option<&String> {
        let children = self.tree.parent(index);
        match children {
            Some(x) => self.tree.get(x),
            None => None,
        }
    }

    pub fn intersect(&self, group_index: Index, ray: &Ray) -> Vec<Intersection> {
        let mut intersections = self
            .children(group_index)
            .iter()
            .map(|x| {
                let shape = self.shapes.get(x.clone()).unwrap();
                shape.intersect_group(ray, shape.transformation())
            })
            .flatten()
            .collect::<Vec<_>>();
        intersections.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
        intersections
    }

    pub fn world_to_object(&self, shape_index: Index, point: Vector4<f32>) -> Vector4<f32> {
        let shape = self.get_shape_by_index(shape_index);
        let inv_op = shape.unwrap().transformation().try_inverse();
        let inv = inv_op.unwrap();
        let parent = self.tree.parent(shape_index);
        match parent {
            Some(_) => inv * self.world_to_object(parent.unwrap(), point),
            None => inv * point,
        }
    }

    pub fn normal_to_world(&self, shape_index: Index, normal: Vector4<f32>) -> Vector4<f32> {
        let shape = self.get_shape_by_index(shape_index).unwrap();
        let mut normal = shape.transformation().try_inverse().unwrap().transpose() * normal;
        normal.w = 0.;
        normal = normal.normalize();
        let parent = self.tree.parent(shape_index);
        match parent {
            Some(parent_index) => self.normal_to_world(parent_index, normal),
            None => normal,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Group {
    pub id: String,
    pub transformation: Matrix4<f32>,
    pub material: Material,
}

impl PartialEq for Group {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

#[allow(dead_code)]
impl Group {
    pub fn new(transformation: Option<Matrix4<f32>>, material: Option<Material>) -> Group {
        Group {
            id: Shape::generate_id(Some("group")),
            transformation: match transformation {
                Some(x) => x,
                None => Matrix4::<f32>::identity(),
            },
            material: match material {
                Some(x) => x,
                None => Material::new_default(),
            },
        }
    }

    pub fn new_default() -> Group {
        Group::new(None, None)
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use nalgebra::Vector4;

    use crate::{
        epsilon::EPSILON, matrix::CGMatrix, plane::Plane, sphere::Sphere, tuple::TupleOperation,
    };

    use std::f32::consts::PI;

    use super::*;
    #[test]
    fn create_new_group() {
        let g = Group::new_default();
        // assert_eq!(g.children.len(), 0);
    }

    #[test]
    fn shape_has_parent_attribute() {
        let s = Shape::Plane(Plane::new_default());
        // assert_eq!(s.parent(), None);
    }

    #[test]
    fn add_child_to_group() {
        let g = Shape::Group(Group::new_default());
        let s = Shape::Plane(Plane::new_default());

        let mut tree = VecTree::new();
        let gn = tree.insert_root(g.id());
        let sn = tree.insert(s.id(), gn);
        let scene = Scene::from_shapes(tree.clone(), vec![g.clone(), s.clone()]);

        assert_eq!(scene.len_children(gn), 1);
        let child = *scene.children(gn).get(0).unwrap();
        assert_eq!(child, &s.id());
        let parent = scene.parent(sn).unwrap();
        assert_eq!(parent, &g.id());
    }

    #[test]
    fn intersectiong_ray_with_nonempty_group() {
        let g = Shape::Group(Group::new_default());
        let s1 = Shape::Sphere(Sphere::new_default());
        let s2 = Shape::Sphere(Sphere::new(
            Some(Matrix4::translation(0., 0., -3.)),
            None,
            None,
            None,
        ));
        let s3 = Shape::Sphere(Sphere::new(
            Some(Matrix4::translation(5., 0., 0.)),
            None,
            None,
            None,
        ));

        let mut tree = VecTree::new();
        let gn = tree.insert_root(g.id());
        tree.insert(s1.id(), gn);
        tree.insert(s2.id(), gn);
        tree.insert(s3.id(), gn);
        let scene = Scene::from_shapes(tree, vec![g.clone(), s1.clone(), s2.clone(), s3.clone()]);

        let ray = Ray::from_tuple((0., 0., -5.), (0., 0., 1.));
        let xs = scene.intersect(gn, &ray);

        assert_eq!(xs.len(), 4);
        assert_eq!(xs[0].shape, &s2);
        assert_eq!(xs[1].shape, &s2);
        assert_eq!(xs[2].shape, &s1);
        assert_eq!(xs[3].shape, &s1);
    }

    #[test]
    fn intersectiong_transformed_group() {
        let g = Shape::Group(Group::new(Some(Matrix4::scaling(2., 2., 2.)), None));
        let s = Shape::Sphere(Sphere::new(
            Some(Matrix4::translation(5., 0., 0.)),
            None,
            None,
            None,
        ));
        let mut tree = VecTree::new();
        let gn = tree.insert_root(g.id());
        tree.insert(s.id(), gn);
        let scene = Scene::from_shapes(tree, vec![g.clone(), s.clone()]);
        let ray = Ray::from_tuple((10., 0., -10.), (0., 0., 1.));
        let xs = scene.intersect(gn, &ray);

        assert_eq!(xs.len(), 2);
    }

    #[test]
    fn converting_point_from_world_to_object_space() {
        let g1 = Shape::Group(Group::new(Some(Matrix4::rotation_y(PI / 2.0)), None));
        let g2 = Shape::Group(Group::new(Some(Matrix4::scaling(2., 2., 2.)), None));
        let s = Shape::Sphere(Sphere::new(
            Some(Matrix4::translation(5., 0., 0.)),
            None,
            None,
            None,
        ));
        let mut tree = VecTree::new();
        let g1n = tree.insert_root(g1.id());
        let g2n = tree.insert(g2.id(), g1n);
        let sn = tree.insert(s.id(), g2n);
        let scene = Scene::from_shapes(tree, vec![g1.clone(), g2.clone(), s.clone()]);
        let p = scene.world_to_object(sn, Vector4::point(-2., 0., -10.));
        assert_relative_eq!(p, Vector4::point(0., 0., -1.), epsilon = EPSILON)
    }

    #[test]
    fn converting_normal_from_object_to_world_space() {
        let g1 = Shape::Group(Group::new(Some(Matrix4::rotation_y(PI / 2.0)), None));
        let g2 = Shape::Group(Group::new(Some(Matrix4::scaling(1., 2., 3.)), None));
        let s = Shape::Sphere(Sphere::new(
            Some(Matrix4::translation(5., 0., 0.)),
            None,
            None,
            None,
        ));
        let mut tree = VecTree::new();
        let g1n = tree.insert_root(g1.id());
        let g2n = tree.insert(g2.id(), g1n);
        let sn = tree.insert(s.id(), g2n);
        let scene = Scene::from_shapes(tree, vec![g1.clone(), g2.clone(), s.clone()]);
        let n = scene.normal_to_world(
            sn,
            Vector4::point(
                3.0_f32.sqrt() / 3.0,
                3.0_f32.sqrt() / 3.0,
                3.0_f32.sqrt() / 3.0,
            ),
        );
        assert_relative_eq!(
            n,
            Vector4::point(0.2857, 0.4286, -0.8571),
            epsilon = EPSILON * 10.
        )
    }
}
