//// A fourth order Runge-Kutta method

fn order1(fun: &dyn Fn(f64,f64,f64)->f64,x:f64,y:f64,z:f64,h:f64) -> f64 {
    h*fun(x,y,z)
}

fn order2(fun: &dyn Fn(f64,f64,f64)->f64,x:f64,y:f64,z:f64,h:f64,k1:f64,l1:f64) -> f64 {
    h*fun(x+0.5*h,y+0.5*k1,z+0.5*l1)
}
fn order3(fun: &dyn Fn(f64,f64,f64)->f64,x:f64,y:f64,z:f64,h:f64,k2:f64,l2:f64) -> f64 {
    h*fun(x+0.5*h,y+0.5*k2,z+0.5*l2)
}
fn order4(fun: &dyn Fn(f64,f64,f64)->f64,x:f64,y:f64,z:f64,h:f64,k3:f64,l3:f64) -> f64 {
    h*fun(x+h,y+k3,z+l3)
}

fn dy_and_dz(
    yprime: &dyn Fn(f64,f64,f64)->f64,
    zprime: &dyn Fn(f64,f64,f64)->f64,
    x:f64, y:f64, z:f64, h:f64
) -> (f64,f64) {
    let k1 = order1(yprime, x, y, z, h);
    let l1 = order1(zprime, x, y, z, h);
    let k2 = order2(yprime, x, y, z, h, k1, l1);
    let l2 = order2(zprime, x, y, z, h, k1, l1);
    let k3 = order3(yprime, x, y, z, h, k2, l2);
    let l3 = order3(zprime, x, y, z, h, k2, l2);
    let k4 = order4(yprime, x, y, z, h, k3, l3);
    let l4 = order4(zprime, x, y, z, h, k3, l3);
    (
        k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0,
        l1/6.0 + l2/3.0 + l3/3.0 + l4/6.0
    )
}

pub fn get_next_xyz(
    yprime: &dyn Fn(f64,f64,f64)->f64,
    zprime: &dyn Fn(f64,f64,f64)->f64,
    x:f64, y:f64, z:f64, h:f64
) -> (f64,f64,f64) {
    let (dy,dz) = dy_and_dz(yprime, zprime, x, y, z, h);
    (x+h, y+dy, z+dz)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn add(x:f64, y:f64, z:f64) -> f64 {
        x+y+z
    }
    #[test]
    fn test_order1() {
        let k1:f64 = order1(&add, 1.0, 2.0, 3.0, 2.0);
        assert_eq!(k1, 12.0, "k1");
    }

}