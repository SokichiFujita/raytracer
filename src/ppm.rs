extern crate approx;
extern crate ulid;

use crate::{canvas::Canvas, color::Color};
use std::fs::File;
use std::io::{Error, Write};

#[derive(Clone, Debug)]
pub struct Ppm {
    header: String,
    body: String,
}

impl Ppm {
    pub fn header(width: usize, height: usize) -> String {
        format!("P3\n{} {}\n255\n", width, height)
    }

    pub fn new_empty(width: usize, height: usize) -> Self {
        Ppm {
            header: Ppm::header(width, height),
            body: String::new(),
        }
    }

    pub fn from_canvas(canvas: Canvas) -> Self {
        Ppm::new(canvas.width, canvas.width, &canvas.pixels)
    }

    pub fn new(width: usize, height: usize, pixels: &Vec<Color>) -> Self {
        const PIXEL_PER_LINE: usize = 4;

        let max_body_line_num = pixels.len() / PIXEL_PER_LINE as usize;
        let mut body = String::new().to_owned();
        for n in 0..max_body_line_num + 1 {
            for m in 0..PIXEL_PER_LINE {
                let a = pixels.get(n * PIXEL_PER_LINE + m);
                match a {
                    Some(_) => {
                        let color = a.unwrap().to_u8();
                        body.push_str(&color.r.to_string());
                        body.push_str(" ");
                        body.push_str(&color.g.to_string());
                        body.push_str(" ");
                        body.push_str(&color.b.to_string());
                        body.push_str(" ");
                    }
                    None => {
                        //Nothing to do
                    }
                }
            }
            body.push_str("\n");
        }
        Ppm {
            header: Ppm::header(width, height),
            body: body,
        }
    }

    pub fn render(&self, filename: &'static str) -> Result<usize, Error> {
        let path = format!("ppms/{}.ppm", filename);
        let mut f = File::create(&path)?;
        f.write_all(self.header.as_bytes())?;
        f.write_all(self.body.as_bytes())?;
        println!("Completed: /ppms/{:}.ppm", filename.to_string());
        Ok(0)
    }

    pub fn out(&self, filename: &'static str) {
        let result = self.render(filename);
        if result.is_err() {
            println!("Error: /ppms/{:?}.ppm", filename);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn ppm_out() {
        let color = Color::rgb(1.0, 0.6, 0.4);
        let color2 = Color::rgb(0.1, 0.2, 0.8);
        let colors = vec![
            color.clone(),
            color2.clone(),
            color.clone(),
            color2.clone(),
            color.clone(),
            color2.clone(),
        ];
        let width = 2;
        let height = 3;
        assert_eq!(colors.len(), width * height);
        let ppm = Ppm::new(width, height, &colors);
        assert_eq!(
            ppm.header.len() > 5 + width.to_string().len() + height.to_string().len(),
            true
        );
        assert_eq!(ppm.body.len() > width * 3 * 4, true);
        ppm.out("canvas_test");
    }
}
