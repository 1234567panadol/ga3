mod ga3;

use ga3::*;

fn main() {
    let vec = OddVersor::from_vec(1.0, 0.0, 0.0);
    let bivec = EvenVersor::from_bivec(0.0, 1.0, 0.0);
    let rotor = Rotor::from_bivector_angle(bivec, 45.0f64.to_radians());
    let result = rotor.sandwich(vec);
    
    println!("{}", result);
}