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

/// 2^e в Q16.16 (точное для целых e в домене |e|<=31).
pub fn q_pow2(e: i32) -> Option<Q> {
    if !(-31..=31).contains(&e) {
        return None;
    }
    if e >= 0 {
        Some(1i64 << (16 + e))
    } else {
        Some((1i64 << 16) >> (-e))
    }
}

/// floor(log2(x)) для x>0 в Q16.16 (домен: x строго положителен).
pub fn q_log2_i(x: Q) -> Option<i32> {
    if x <= 0 {
        return None;
    }
    let v = x as u64;
    Some((63 - v.leading_zeros()) as i32 - 16)
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

    #[test]
    fn pow2_and_log2() {
        assert_eq!(q_pow2(0), Some(q_from_i32(1)));
        assert_eq!(q_pow2(1), Some(q_from_i32(2)));
        assert_eq!(q_pow2(-1), Some(32768));
        assert_eq!(q_log2_i(q_from_i32(8)), Some(3));
        assert_eq!(q_log2_i(q_from_i32(1)), Some(0));
        assert!(q_log2_i(0).is_none());
        assert!(q_pow2(100).is_none());
        for e in -8..=8 {
            let p = q_pow2(e).unwrap();
            assert_eq!(q_log2_i(p), Some(e), "e={e}");
        }
    }
}
