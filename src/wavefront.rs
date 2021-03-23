extern crate approx;
extern crate ulid;

use crate::{matrix::CGMatrix, tuple::TupleOperation};
use na::{Matrix4, Vector4};
use std::fs::File;
use std::io::{BufRead, BufReader, Error, Read, Write};

#[derive(Clone, Debug)]
pub struct Wavefront {
    comments: Vec<String>,
    vertices: Vec<Vector4<f32>>,
    faces: Vec<Vector4<f32>>,
    groups: Vec<Vector4<f32>>,
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
        }
    }

    pub fn parse(&mut self, filename: &'static str) -> Result<usize, Error> {
        let path = format!("wavefront_obj_files/{}", filename);
        for result in BufReader::new(File::open(path)?).lines() {
            let line = result?;
            let words: Vec<&str> = line.split_whitespace().collect();
            if words.len() == 0 {
                self.comments.push(line);
                break;
            }
            if words[0].len() != 1 {
                self.comments.push(line);
                break;
            }
            let category = words[0].chars().nth(0).unwrap();
            let tuple = Vector4::point(
                words[1].parse::<f32>().unwrap(),
                words[2].parse::<f32>().unwrap(),
                words[3].parse::<f32>().unwrap(),
            );
            match category {
                'v' => self.vertices.push(tuple),
                'f' => self.faces.push(tuple),
                'g' => self.groups.push(tuple),
                _ => self.comments.push(line),
            }
        }
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
    fn file_read() {
        let mut w = Wavefront::new_empty();
        w.parse("vertex.obj");
        assert_eq!(w.vertices[0], Vector4::point(-1., 1., 0.));
        assert_eq!(w.vertices[1], Vector4::point(-1., 0.5, 0.));
        assert_eq!(w.vertices[2], Vector4::point(1., 0., 0.));
        assert_eq!(w.vertices[3], Vector4::point(1., 1., 0.));
    }
}
