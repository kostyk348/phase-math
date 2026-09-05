//! Q16.16 фиксированная точка с домен-валидацией на границе фазы.
//! Операции возвращают Option: нарушение домена -> None -> вызывающий
//! откатывает фазу (fail = un-run).

pub type Q = i64; // Q16.16: значение = целое / 65536

pub const FRAC: i64 = 65536;
pub const Q_MAX: Q = i64::MAX;

pub fn q_from_i32(v: i32) -> Q {
    (v as i64) << 16
}

pub fn q_mul(a: Q, b: Q) -> Option<Q> {
    let p = (a as i128) * (b as i128);
    if p < (i64::MIN as i128) || p > (i64::MAX as i128) {
        return None;
    }
    let r = p >> 16;
    Some(r as i64)
}

pub fn q_div(a: Q, b: Q) -> Option<Q> {
    if b == 0 {
        return None;
    }
    let p = ((a as i128) << 16) / (b as i128);
    if p < (i64::MIN as i128) || p > (i64::MAX as i128) {
        return None;
    }
    Some(p as i64)
}

/// Целочисленный sqrt в Q16.16: result_Q = isqrt(a)*256 (домен a>=0).
pub fn q_sqrt(a: Q) -> Option<Q> {
    if a < 0 {
        return None;
    }
    Some((isqrt(a as u64) << 8) as Q)
}

/// Целочисленный квадратный корень (Newton, 64-бит).
fn isqrt(n: u64) -> u64 {
    if n < 2 {
        return n;
    }
    let mut x = n;
    let mut y = (x + 1) / 2;
    while y < x {
        x = y;
        y = (x + n / x) / 2;
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sqrt_exact_squares() {
        for v in [0u64, 1, 4, 9, 16, 25, 100, 65536, 1 << 30] {
            let a = (v << 16) as Q;
            let s = q_sqrt(a).unwrap();
            // (s/65536)^2 ~ v: проверяем через q_mul
            let sq = q_mul(s, s).unwrap();
            let err = (sq - a).abs();
            assert!(err < (1 << 18), "v={v} err={err}");
        }
    }

    #[test]
    fn sqrt_domain() {
        assert!(q_sqrt(-1).is_none());
    }
}
