//! Legacy Q16.16 fixed-point helpers with boundary domain validation.
//!
//! `Q` is a raw `i64`; the represented value is `raw / 2^16`. Operations
//! return `Option`: a domain violation yields `None` and the caller rolls
//! the phase back (fail = un-run). These free functions predate the generic
//! [`crate::Fixed`] API and are kept for compatibility; see
//! `docs/math-spec.md` for the full per-operation error budgets and the
//! (deliberate) rounding differences from `Fixed` (mul floors here, but
//! `q_div` truncates toward zero while `Fixed::checked_div` floors).

/// Q16.16 raw value: real value = `raw / 65536`.
pub type Q = i64;

/// The Q16.16 fractional scale: 2^16 = 65536 subunits per unit
/// (value = `raw / FRAC`).
pub const FRAC: i64 = 65536;
/// The largest representable Q16.16 raw value (equal to `i64::MAX`).
pub const Q_MAX: Q = i64::MAX;

/// Converts an integer to Q16.16 (raw `v << 16`). Exact, injective; the
/// exact inverse on its image is `raw >> 16`. Always succeeds for `i32`.
pub fn q_from_i32(v: i32) -> Q {
    (v as i64) << 16
}

/// Q16.16 multiplication with floor rounding.
///
/// Computes `floor(a*b / 2^16)` in i128. Error bound (exact, worst case):
/// `-2^-16 < err <= 0` (at most 1 ulp, toward -infinity; 0 when the raw
/// product is divisible by `2^16`). Domain (conservative): returns `None`
/// when the raw product `a*b` itself does not fit in `i64` — strictly more
/// restrictive than result representability. Class: lossy.
pub fn q_mul(a: Q, b: Q) -> Option<Q> {
    let p = (a as i128) * (b as i128);
    if p < (i64::MIN as i128) || p > (i64::MAX as i128) {
        return None;
    }
    let r = p >> 16;
    Some(r as i64)
}

/// Q16.16 division with truncation toward zero.
///
/// Computes `(a << 16) / b` in i128. Error bound (exact, worst case):
/// `|err| < 2^-16` (at most 1 ulp; 0 when the division is exact).
/// Domain: `b != 0` and the shifted quotient must fit in `i64`.
/// Class: lossy. Note: `Fixed::checked_div` floors instead — documented
/// deliberately in `docs/math-spec.md`.
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

/// Integer square root in Q16.16 with floor rounding: `isqrt(a) << 8`.
///
/// Error bound (exact, worst case): `-2^-8 < err <= 0` (at most 256 ulp —
/// squaring halves the precision). Domain: `a >= 0`. Class: lossy.
pub fn q_sqrt(a: Q) -> Option<Q> {
    if a < 0 {
        return None;
    }
    Some((isqrt_u64(a as u64) << 8) as Q)
}

/// `2^e` in Q16.16.
///
/// Domain: `e` in `[-31, 31]`, else `None`. Exact (error 0) for
/// `e` in `[-16, 31]`; for `e` in `[-31, -17]` the exact raw `2^(16+e)` is
/// below 1 ulp and the result is `0` with error `2^e < 2^-16`.
/// Class: partially invertible — `q_log2_i` is the exact inverse exactly on
/// `e` in `[-16, 31]`.
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

/// `floor(log2(x))` for Q16.16.
///
/// Returns the exact floor of the real `log2` (error `-1 < err_log2 <= 0`
/// in bits). Domain: `x > 0`, else `None`. Class: lossy (many-to-one); the
/// exact inverse of `q_pow2` on `e` in `[-16, 31]`.
pub fn q_log2_i(x: Q) -> Option<i32> {
    if x <= 0 {
        return None;
    }
    let v = x as u64;
    Some((63 - v.leading_zeros()) as i32 - 16)
}

use crate::fixed::isqrt_u64;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sqrt_exact_squares() {
        for v in [0u64, 1, 4, 9, 16, 25, 100, 65536, 1 << 30] {
            let a = (v << 16) as Q;
            let s = q_sqrt(a).unwrap();
            // (s/65536)^2 ~ v: check via q_mul
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
