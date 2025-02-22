use std::{
    fmt::Display,
    ops::{Add, Mul, Sub},
};

pub trait Versor: Sized + Copy + Add + Sub + Mul<Self::Even> + Mul<Self::Odd> {
    type Even: Versor;
    type Odd: Versor;
    fn mag2(self) -> f64;
    fn rev(self) -> Self;
    fn inv(self) -> Self;
    fn normalize(self) -> Self;

    // Add the sandwich operation
    fn sandwich<V, U>(self, versor: V) -> V
    where
        Self: Mul<V, Output = U>,
        U: Mul<Self, Output = V>,
    {
        self * versor * self.inv()
    }
}

#[derive(Clone, Copy)]
pub struct EvenVersor {
    pub sca: f64,
    pub yz: f64,
    pub zx: f64,
    pub xy: f64,
}

#[derive(Clone, Copy)]
pub struct OddVersor {
    pub xyz: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Versor for EvenVersor {
    type Even = EvenVersor;
    type Odd = OddVersor;
    fn mag2(self) -> f64 {
        self.sca * self.sca + self.yz * self.yz + self.zx * self.zx + self.xy * self.xy
    }
    fn rev(self) -> EvenVersor {
        EvenVersor {
            sca: self.sca,
            yz: -self.yz,
            zx: -self.zx,
            xy: -self.xy,
        }
    }
    fn inv(self) -> EvenVersor {
        self.rev() * EvenVersor::from_sca(1.0 / self.mag2())
    }
    fn normalize(self) -> EvenVersor {
        self * EvenVersor::from_sca(1.0 / self.mag2().sqrt())
    }
}

impl Add for EvenVersor {
    type Output = EvenVersor;
    fn add(self, rhs: Self) -> Self::Output {
        EvenVersor {
            sca: self.sca + rhs.sca,
            yz: self.yz + rhs.yz,
            zx: self.zx + rhs.zx,
            xy: self.xy + rhs.xy,
        }
    }
}

impl Sub for EvenVersor {
    type Output = EvenVersor;
    fn sub(self, rhs: Self) -> Self::Output {
        EvenVersor {
            sca: self.sca - rhs.sca,
            yz: self.yz - rhs.yz,
            zx: self.zx - rhs.zx,
            xy: self.xy - rhs.xy,
        }
    }
}

impl Mul<EvenVersor> for EvenVersor {
    type Output = EvenVersor;
    fn mul(self, rhs: EvenVersor) -> Self::Output {
        EvenVersor {
            sca: self.sca * rhs.sca - self.yz * rhs.yz - self.zx * rhs.zx - self.xy * rhs.xy,
            yz: self.yz * rhs.sca + self.sca * rhs.yz + self.xy * rhs.zx - self.zx * rhs.xy,
            zx: self.zx * rhs.sca + self.sca * rhs.zx + self.yz * rhs.xy - self.xy * rhs.yz,
            xy: self.xy * rhs.sca + self.sca * rhs.xy + self.zx * rhs.yz - self.yz * rhs.zx,
        }
    }
}

impl Mul<OddVersor> for EvenVersor {
    type Output = OddVersor;
    fn mul(self, rhs: OddVersor) -> Self::Output {
        OddVersor {
            xyz: self.sca * rhs.xyz + self.yz * rhs.x + self.zx * rhs.y + self.xy * rhs.z,
            x: self.sca * rhs.x - self.yz * rhs.xyz + self.xy * rhs.y - self.zx * rhs.z,
            y: self.sca * rhs.y - self.zx * rhs.xyz + self.yz * rhs.z - self.xy * rhs.x,
            z: self.sca * rhs.z - self.xy * rhs.xyz + self.zx * rhs.x - self.yz * rhs.y,
        }
    }
}

impl EvenVersor {
    pub fn from_sca(sca: f64) -> EvenVersor {
        EvenVersor {
            sca,
            yz: 0.0,
            zx: 0.0,
            xy: 0.0,
        }
    }

    pub fn from_bivec(yz: f64, zx: f64, xy: f64) -> EvenVersor {
        EvenVersor {
            sca: 0.0,
            yz,
            zx,
            xy,
        }
    }

    pub fn sca(self) -> EvenVersor {
        EvenVersor {
            sca: self.sca,
            yz: 0.0,
            zx: 0.0,
            xy: 0.0,
        }
    }

    pub fn bivec(self) -> EvenVersor {
        EvenVersor {
            sca: 0.0,
            yz: self.yz,
            zx: self.zx,
            xy: self.xy,
        }
    }
}

impl Display for EvenVersor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "sca: {}, yz: {}, zx: {}, xy: {}",
            self.sca, self.yz, self.zx, self.xy,
        )
    }
}

impl Versor for OddVersor {
    type Even = EvenVersor;
    type Odd = OddVersor;
    fn mag2(self) -> f64 {
        self.xyz * self.xyz + self.x * self.x + self.y * self.y + self.z * self.z
    }
    fn rev(self) -> OddVersor {
        OddVersor {
            xyz: -self.xyz,
            x: self.x,
            y: self.y,
            z: self.z,
        }
    }
    fn inv(self) -> OddVersor {
        self.rev() * EvenVersor::from_sca(1.0 / self.mag2())
    }
    fn normalize(self) -> OddVersor {
        self * EvenVersor::from_sca(1.0 / self.mag2().sqrt())
    }
}

impl Add for OddVersor {
    type Output = OddVersor;
    fn add(self, rhs: Self) -> Self::Output {
        OddVersor {
            xyz: self.xyz + rhs.xyz,
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl Sub for OddVersor {
    type Output = OddVersor;
    fn sub(self, rhs: Self) -> Self::Output {
        OddVersor {
            xyz: self.xyz - rhs.xyz,
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl Mul<OddVersor> for OddVersor {
    type Output = EvenVersor;
    fn mul(self, rhs: OddVersor) -> Self::Output {
        EvenVersor {
            sca: self.x * rhs.x + self.y * rhs.y + self.z * rhs.z - self.xyz * rhs.xyz,
            yz: self.x * rhs.xyz + self.xyz * rhs.x + self.y * rhs.z - self.z * rhs.y,
            zx: self.y * rhs.xyz + self.xyz * rhs.y + self.z * rhs.x - self.x * rhs.z,
            xy: self.z * rhs.xyz + self.xyz * rhs.z + self.x * rhs.y - self.y * rhs.x,
        }
    }
}

impl Mul<EvenVersor> for OddVersor {
    type Output = OddVersor;
    fn mul(self, rhs: EvenVersor) -> Self::Output {
        OddVersor {
            xyz: self.xyz * rhs.sca + self.x * rhs.yz + self.y * rhs.zx + self.z * rhs.xy,
            x: self.x * rhs.sca - self.xyz * rhs.yz + self.z * rhs.zx - self.y * rhs.xy,
            y: self.y * rhs.sca - self.xyz * rhs.zx + self.x * rhs.xy - self.z * rhs.yz,
            z: self.z * rhs.sca - self.xyz * rhs.xy + self.y * rhs.yz - self.x * rhs.zx,
        }
    }
}

impl OddVersor {
    pub fn from_vec(x: f64, y: f64, z: f64) -> OddVersor {
        OddVersor { xyz: 0.0, x, y, z }
    }

    pub fn vec(self) -> OddVersor {
        OddVersor {
            xyz: 0.0,
            x: self.x,
            y: self.y,
            z: self.z,
        }
    }
}

impl Display for OddVersor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "xyz: {}, x: {}, y: {}, z: {}",
            self.xyz, self.x, self.y, self.z,
        )
    }
}

pub type Rotor = EvenVersor;

impl Rotor {
    pub fn from_vectors(a: OddVersor, b: OddVersor) -> Rotor {
        let vec_a = a.vec();
        let vec_b = b.vec();
        let normalized_a = vec_a.normalize();
        let normalized_b = vec_b.normalize();
        let mid = (normalized_a + normalized_b).normalize();
        normalized_a * mid
    }
    
    pub fn from_bivector_angle(even_versor: EvenVersor, angle: f64) -> Rotor {
        let bivec = even_versor.bivec();
        let angle_half = angle / 2.0;
        EvenVersor::from_sca(angle_half.cos()) + bivec.normalize() * EvenVersor::from_sca(angle_half.sin())
    }
}
