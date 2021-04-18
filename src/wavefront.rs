extern crate approx;
extern crate ulid;

use crate::{
    group::Group, group::Scene, matrix::CGMatrix, shape::Shape, smooth_triangle::SmoothTriangle,
    triangle::Triangle, tuple::TupleOperation,
};
use na::{Matrix4, Vector4};
use std::fs::File;
use std::io::{BufRead, BufReader, Error, Read, Write};
use vec_tree::{Index, VecTree};

#[derive(Clone, Debug)]
pub struct Wavefront {
    pub comments: Vec<String>,
    pub vertices: Vec<Vector4<f32>>,
    pub faces: Vec<Vector4<f32>>,
    pub groups: Vec<String>,
    pub scene: Option<Scene>,
    pub normals: Vec<Vector4<f32>>,
}

impl Wavefront {
    pub fn header(width: usize, height: usize) -> String {
        format!("P3\n{} {}\n255\n", width, height)
    }

    pub fn new_empty() -> Self {
        Wavefront {
            comments: vec![],
            vertices: vec![],
            faces: vec![],
            groups: vec![],
            scene: None,
            normals: vec![],
        }
    }

    pub fn parse(&mut self, filename: &'static str) -> Result<usize, Error> {
        let path = format!("wavefront_obj_files/{}", filename);
        let g = Shape::Group(Group::new_default());

        let mut tree = VecTree::new();
        let root = tree.insert_root(g.id());

        let mut shapes = vec![];

        self.vertices.push(Vector4::point(0., 0., 0.));
        self.normals.push(Vector4::point(0., 0., 0.));

        for result in BufReader::new(File::open(path)?).lines() {
            let line = result?;

            let mut words: Vec<&str> = line.split_whitespace().collect();
            if words.len() == 0 {
                self.comments.push(line);
                continue;
            }
            if words[0].len() > 2 {
                self.comments.push(line);
                continue;
            }
            let category = words[0]; //first word
            match category {
                "v" => {
                    let tuple = Vector4::point(
                        words[1].parse::<f32>().unwrap(),
                        words[2].parse::<f32>().unwrap(),
                        words[3].parse::<f32>().unwrap(),
                    );
                    self.vertices.push(tuple);
                }
                "g" => {
                    self.groups.push(line);
                }
                "vn" => {
                    let tuple = Vector4::vector(
                        words[1].parse::<f32>().unwrap(),
                        words[2].parse::<f32>().unwrap(),
                        words[3].parse::<f32>().unwrap(),
                    );
                    self.normals.push(tuple);
                }
                "f" => {
                    if words[1].contains("/") {
                        words.remove(0);
                        let mut vertices = vec![];
                        let mut normals = vec![];

                        for word in words {
                            let f = word.split("/").collect::<Vec<&str>>();
                            let ii1 = f[0].parse::<usize>().unwrap();
                            let ii3 = f[2].parse::<usize>().unwrap();
                            let p = self.vertices[ii1];
                            let n = self.normals[ii3];
                            vertices.push(p);
                            normals.push(n)
                        }

                        let l = vertices.len();
                        if l > 3 {
                            //Create multiple triangle polygons from N(>3) vertices
                            for index in 2..l {
                                let p1 = self.vertices[1];
                                let p2 = self.vertices[index];
                                let p3 = self.vertices[index + 1];
                                let n1 = self.normals[1];
                                let n2 = self.normals[index];
                                let n3 = self.normals[index + 1];

                                let triangle = Shape::SmoothTriangle(SmoothTriangle::new_with_n(
                                    None, None, p1, p2, p3, n1, n2, n3,
                                ));
                                tree.insert(triangle.id(), root);
                                shapes.push(triangle);
                            }
                        } else if l == 3 {
                            let triangle = Shape::SmoothTriangle(SmoothTriangle::new_with_n(
                                None,
                                None,
                                vertices[0],
                                vertices[1],
                                vertices[2],
                                normals[0],
                                normals[1],
                                normals[2],
                            ));
                            tree.insert(triangle.id(), root);
                            shapes.push(triangle);
                        } else {
                            //Error
                        }
                    } else {
                        // Create one triangle polygon from 3 vertices
                        words.remove(0);
                        let mut vertices = vec![];
                        for word in words {
                            let ii = word.parse::<usize>().unwrap();
                            let pi = self.vertices[ii];
                            vertices.push(pi);
                        }
                        let l = vertices.len();
                        if l > 3 {
                            //Create multiple triangle polygons from N(>3) vertices
                            for index in 2..l {
                                let p1 = self.vertices[1];
                                let p2 = self.vertices[index];
                                let p3 = self.vertices[index + 1];

                                let triangle = Shape::SmoothTriangle(SmoothTriangle::new(
                                    None, None, p1, p2, p3,
                                ));
                                tree.insert(triangle.id(), root);
                                shapes.push(triangle);
                            }
                        } else if l == 3 {
                            // Create one triangle polygon from 3 vertices
                            let triangle = Shape::SmoothTriangle(SmoothTriangle::new(
                                None,
                                None,
                                vertices[0],
                                vertices[1],
                                vertices[2],
                            ));
                            tree.insert(triangle.id(), root);
                            shapes.push(triangle);
                        } else {
                            // Error
                        }
                    };
                }
                _ => self.comments.push(line),
            }
        }

        let scene = Scene::from_shapes(tree, shapes);
        self.scene = Some(scene);

        println!(
            "Completed: /wavefront_obj_files/{:}.obj",
            filename.to_string()
        );
        Ok(0)
    }

    pub fn out(&self, filename: &'static str) {}
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn parse_vertex() {
        let mut w = Wavefront::new_empty();
        let result = w.parse("vertex.obj");
        assert_eq!(result.is_ok(), true);
        assert_eq!(w.vertices[1], Vector4::point(-1., 1., 0.));
        assert_eq!(w.vertices[2], Vector4::point(-1., 0.5, 0.));
        assert_eq!(w.vertices[3], Vector4::point(1., 0., 0.));
        assert_eq!(w.vertices[4], Vector4::point(1., 1., 0.));
    }

    #[test]
    fn parse_face() {
        let mut w = Wavefront::new_empty();
        let result = w.parse("face.obj");
        assert_eq!(result.is_ok(), true);

        let t = w.scene.unwrap();
        let shapes = t.shapes;
        let tree = t.tree;
        let root = tree.get_root_index().unwrap();

        let nodes: Vec<Index> = tree.children(root).collect();
        let i1 = tree.get(nodes[0]).unwrap();
        let s1 = shapes.get(i1).unwrap();

        assert_eq!(s1.smooth_triangle().unwrap().p1, w.vertices[1]);
        assert_eq!(s1.smooth_triangle().unwrap().p2, w.vertices[2]);
        assert_eq!(s1.smooth_triangle().unwrap().p3, w.vertices[3]);

        let i2 = tree.get(nodes[1]).unwrap();
        let s2 = shapes.get(i2).unwrap();

        assert_eq!(s2.smooth_triangle().unwrap().p1, w.vertices[1]);
        assert_eq!(s2.smooth_triangle().unwrap().p2, w.vertices[3]);
        assert_eq!(s2.smooth_triangle().unwrap().p3, w.vertices[4]);
    }

    #[test]
    fn parse_polygone() {
        let mut w = Wavefront::new_empty();
        let result = w.parse("polygones.obj");
        assert_eq!(result.is_ok(), true);

        let t = w.scene.unwrap();
        let shapes = t.shapes;
        let tree = t.tree;
        let root = tree.get_root_index().unwrap();

        let nodes: Vec<Index> = tree.children(root).collect();
        let i1 = tree.get(nodes[0]).unwrap();
        let s1 = shapes.get(i1).unwrap();

        assert_eq!(s1.smooth_triangle().unwrap().p1, w.vertices[1]);
        assert_eq!(s1.smooth_triangle().unwrap().p2, w.vertices[2]);
        assert_eq!(s1.smooth_triangle().unwrap().p3, w.vertices[3]);

        let i2 = tree.get(nodes[1]).unwrap();
        let s2 = shapes.get(i2).unwrap();

        assert_eq!(s2.smooth_triangle().unwrap().p1, w.vertices[1]);
        assert_eq!(s2.smooth_triangle().unwrap().p2, w.vertices[3]);
        assert_eq!(s2.smooth_triangle().unwrap().p3, w.vertices[4]);

        let i3 = tree.get(nodes[2]).unwrap();
        let s3 = shapes.get(i3).unwrap();

        assert_eq!(s3.smooth_triangle().unwrap().p1, w.vertices[1]);
        assert_eq!(s3.smooth_triangle().unwrap().p2, w.vertices[4]);
        assert_eq!(s3.smooth_triangle().unwrap().p3, w.vertices[5]);
    }

    #[test]
    fn parse_vertex_normal() {
        let mut w = Wavefront::new_empty();
        let result = w.parse("vertex_normal.obj");
        assert_eq!(result.is_ok(), true);

        assert_eq!(w.normals[1], Vector4::vector(0., 0., 1.));
        assert_eq!(w.normals[2], Vector4::vector(0.707, 0., -0.707));
        assert_eq!(w.normals[3], Vector4::vector(1., 2., 3.));
    }

    #[test]
    fn parse_face_with_normal() {
        let mut w = Wavefront::new_empty();
        let result = w.parse("face_with_normal.obj");
        assert_eq!(result.is_ok(), true);

        let t = w.scene.unwrap();
        let shapes = t.shapes;
        let tree = t.tree;
        let root = tree.get_root_index().unwrap();

        let nodes: Vec<Index> = tree.children(root).collect();
        let i1 = tree.get(nodes[0]).unwrap();
        let s1 = shapes.get(i1).unwrap();
        let i2 = tree.get(nodes[1]).unwrap();
        let s2 = shapes.get(i2).unwrap();

        assert_eq!(s1.smooth_triangle().unwrap().p1, w.vertices[1]);
        assert_eq!(s1.smooth_triangle().unwrap().p2, w.vertices[2]);
        assert_eq!(s1.smooth_triangle().unwrap().p3, w.vertices[3]);
        assert_eq!(s1.smooth_triangle().unwrap().n1, w.normals[3]);
        assert_eq!(s1.smooth_triangle().unwrap().n2, w.normals[1]);
        assert_eq!(s1.smooth_triangle().unwrap().n3, w.normals[2]);

        assert_eq!(s2.smooth_triangle().unwrap().p1, w.vertices[1]);
        assert_eq!(s2.smooth_triangle().unwrap().p2, w.vertices[2]);
        assert_eq!(s2.smooth_triangle().unwrap().p3, w.vertices[3]);
        assert_eq!(s2.smooth_triangle().unwrap().n1, w.normals[3]);
        assert_eq!(s2.smooth_triangle().unwrap().n2, w.normals[1]);
        assert_eq!(s2.smooth_triangle().unwrap().n3, w.normals[2]);
    }
}
