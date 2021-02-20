pub struct RefractiveIndex;

#[allow(dead_code)]
impl RefractiveIndex {
    pub const VACCUME: f32 = 1.0;
    pub const AIR: f32 = 1.00029;
    pub const WATER: f32 = 1.333;
    pub const GLASS: f32 = 1.5;
    pub const DIAMOND: f32 = 2.417;
}

#[cfg(test)]
mod tests {}
