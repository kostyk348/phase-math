//! # phase-math — reversible math for the PHASE paradigm
//!
//! Math **inside a phase** = bijections on Z/2^64 with analytically computed
//! inverses (composition of any ops = bijection → any formula is exactly
//! undoable: `F⁻¹(F(x)) == x`). Float math is not invertible and stays **at
//! the phase boundary**.
//!
//! Families provided by this crate:
//!
//! - `MulOdd(m)`, `Add`, `Xor`, `RotL` — bijective primitives with
//!   `.apply()` / `.inverse()` (exact, error 0);
//! - `arx_round` / `arx_round_inv` — reversible ARX round on a word pair;
//! - `Compose` — builder for sequences with automatic analytic inversion;
//! - [`Fixed`]`<FRAC>` — fixed-point with an explicit Q-format scale,
//!   checked (returns `Option`) and saturating (returns a value plus a
//!   saturation flag) variants, documented worst-case error bounds;
//! - [`Interval`]`<FRAC>` — closed intervals with outward rounding: the
//!   exact real result is guaranteed to stay inside;
//! - legacy `q_*` free functions over `Q = i64` (Q16.16) with boundary
//!   domain validation (`None` ⇒ the caller rolls the phase back).
//!
//! Per-operation classification (bijective / partially invertible / lossy)
//! and the full error-budget tables are specified — and mapped to tests —
//! in `docs/math-spec.md`.
//!
//! # Examples
//!
//! ```
//! use phase_math::{Compose, Interval, Op, Q16};
//!
//! // 1. Bijections on Z/2^64: any composition inverts analytically.
//! let mut enc = Compose::new();
//! enc.push(Op::MulOdd(0x9E37_79B9_7F4A_7C15));
//! enc.push(Op::Add(0xD1B5_4A32_D192_ED03));
//! assert_eq!(enc.inverse().apply(enc.apply(42)), 42);
//!
//! // 2. Fixed point (Q16.16 here): exact when representable, None when not.
//! let x = Q16::from_i32(3).unwrap();
//! let y = Q16::from_i32(7).unwrap();
//! assert_eq!(x.checked_mul(y), Q16::from_i32(21));
//! assert_eq!(x.checked_mul(Q16::MAX), None); // out of domain -> roll back
//!
//! // 3. Saturating variants report clamping with a flag.
//! let s = Q16::MAX.saturating_add(Q16::ONE);
//! assert!(s.saturated && s.value == Q16::MAX);
//!
//! // 4. Intervals with outward rounding: the exact product stays inside.
//! let around_one = Interval::new(Q16::from_raw(65535), Q16::from_raw(65537)).unwrap();
//! let p = around_one.mul(Interval::point(x)).unwrap(); // 3.0 * (1 ± 1 ulp)
//! assert!(p.contains(x));
//! println!("3.0 * [1 ∓ 1 ulp] = [{}, {}]", p.lo, p.hi);
//! ```

mod biject;
mod fixed;
mod interval;
mod q;

pub use biject::*;
pub use fixed::{Fixed, Q16, Saturated};
pub use interval::Interval;
pub use q::{q_div, q_from_i32, q_log2_i, q_mul, q_pow2, q_sqrt, FRAC, Q, Q_MAX};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compose_roundtrip() {
        let mut c = Compose::new();
        c.push(Op::MulOdd(0x9E3779B97F4A7C15));
        c.push(Op::Add(0xD1B54A32D192ED03));
        c.push(Op::Xor(0xABCDEF0123456789));
        c.push(Op::RotL(19));
        c.push(Op::MulOdd(3));
        let inv = c.inverse();
        let mut rng = 0xF00Du64;
        let mut next = move || {
            rng ^= rng >> 12;
            rng ^= rng << 25;
            rng ^= rng >> 27;
            rng.wrapping_mul(0x2545_F491_4F6C_DD1D)
        };
        for _ in 0..3000 {
            let x = next();
            assert_eq!(inv.apply(c.apply(x)), x);
        }
    }

    #[test]
    fn q_roundtrip_and_domain() {
        // mul: x*y is invertible through division only under preconditions;
        // Q16.16 mul is lossy (floor rounding), but domain-validated.
        // Here: exactness and boundaries.
        let a = q_from_i32(3); // 3.0
        let b = q_from_i32(7); // 7.0
        let p = q_mul(a, b).expect("in domain");
        assert_eq!(p, q_from_i32(21));
        let d = q_div(p, a).expect("ok");
        assert_eq!(d, q_from_i32(7));
        // out of domain -> None (boundary: caller rolls the phase back)
        assert!(q_mul(Q::MAX, Q::MAX).is_none());
        assert!(q_div(a, 0i64.into()).is_none());
        // sqrt exactness for perfect squares
        let s = q_sqrt(q_from_i32(16)).expect("16 -> 4");
        assert_eq!(s, q_from_i32(4));
    }
}
