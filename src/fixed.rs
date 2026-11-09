//! Generic fixed-point arithmetic with an explicit Q-format scale.
//!
//! A [`Fixed`]`<FRAC>` stores a signed 64-bit raw integer; the represented
//! real value is always `raw / 2^FRAC`. The crate is integer-only: there are
//! no floats anywhere, every operation is exactly defined, and every
//! operation carries a worst-case error bound that is derived in
//! `docs/math-spec.md` and verified by randomized tests against exact i128
//! references (see the `tests/` directory).
//!
//! Constraints on `FRAC`:
//!
//! - `FRAC <= 63` (required by division and formatting; larger values
//!   panic in debug builds and wrap in release builds);
//! - even `FRAC` additionally enables [`Fixed::sqrt`];
//! - [`Fixed::ONE`] requires `FRAC <= 62`.
//!
//! The classic instance is [`Q16`] (`FRAC = 16`, historically written
//! "Q16.16": a 64-bit word with 16 fractional bits, i.e. 48 integer bits
//! including the sign).
//!
//! Rounding conventions (deliberate and documented per operation):
//!
//! - [`Fixed::checked_mul`] and [`Fixed::sqrt`] round **down** (floor):
//!   one-sided error, `-2^-FRAC < err <= 0`;
//! - [`Fixed::checked_div`] rounds **down** (floor). This differs from the
//!   legacy `q_div` helper, which truncates toward zero; both behaviors are
//!   documented separately in `docs/math-spec.md`.
//!
//! Saturating variants return a [`Saturated`] value: the (possibly clamped)
//! result plus a flag that is `true` iff clamping actually happened. When
//! the flag is `false` the result is exact; when it is `true` no numeric
//! error bound is promised — the flag is the signal to roll the phase back.

use core::fmt;

/// Result of a saturating operation.
///
/// `value` is the (possibly clamped) result; `saturated` is `true` iff the
/// exact result did not fit and `value` was clamped to `Fixed::MAX` or
/// `Fixed::MIN`. Contract: `saturated == false` implies the result is exact
/// (error 0); `saturated == true` implies `value` is the range clamp and no
/// numeric error bound is offered.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Saturated<T> {
    /// The (possibly clamped) result of the operation.
    pub value: T,
    /// `true` iff the exact result was out of range and `value` was clamped.
    pub saturated: bool,
}

impl<T> Saturated<T> {
    /// Drops the flag, returning the value.
    pub fn into_value(self) -> T {
        self.value
    }
}

/// Signed 64-bit fixed-point number: value = `raw / 2^FRAC`.
///
/// Ordering, equality and hashing follow the raw integer, which is the same
/// as ordering the represented real values (the map is strictly monotone).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct Fixed<const FRAC: u32>(i64);

impl<const FRAC: u32> Fixed<FRAC> {
    /// The value 0 (raw `0`).
    pub const ZERO: Self = Self(0);
    /// The value 1 (raw `2^FRAC`). Requires `FRAC <= 62`.
    pub const ONE: Self = Self(1i64 << FRAC);
    /// The largest representable value (raw `i64::MAX`).
    pub const MAX: Self = Self(i64::MAX);
    /// The smallest representable value (raw `i64::MIN`).
    pub const MIN: Self = Self(i64::MIN);

    /// Wraps a raw integer into a fixed-point number.
    ///
    /// The value of the result is `raw / 2^FRAC`. Infallible: every `i64`
    /// raw is a valid `Fixed`.
    ///
    /// ```
    /// use phase_math::Q16;
    /// assert_eq!(Q16::from_raw(131072), Q16::from_i32(2).unwrap());
    /// ```
    pub const fn from_raw(raw: i64) -> Self {
        Self(raw)
    }

    /// The raw integer behind this value; the real value is `raw() / 2^FRAC`.
    pub const fn raw(self) -> i64 {
        self.0
    }

    /// Converts an integer `v` to `v / 1` in fixed-point (raw `v << FRAC`).
    ///
    /// Returns `None` if `v * 2^FRAC` does not fit in `i64` (possible when
    /// `FRAC > 32`). Exact: error 0. Injective; the exact inverse for
    /// representable values is `raw() >> FRAC` (itself lossy unless the low
    /// `FRAC` raw bits are zero).
    ///
    /// ```
    /// use phase_math::Q16;
    /// assert_eq!(Q16::from_i32(-2), Some(Q16::from_raw(-131072)));
    /// ```
    pub fn from_i32(v: i32) -> Option<Self> {
        debug_assert!(FRAC <= 63);
        let r = (v as i128) << FRAC;
        if r < i64::MIN as i128 || r > i64::MAX as i128 {
            return None;
        }
        Some(Self(r as i64))
    }

    /// Checked addition: exact, `None` on i64 raw overflow.
    ///
    /// Class: partially invertible — the exact inverse of `a.checked_add(b)`
    /// is `a.checked_sub(b)` (error 0), valid whenever the operation
    /// returned `Some` (documented precondition: no raw overflow).
    ///
    /// ```
    /// use phase_math::Q16;
    /// let x = Q16::from_i32(2).unwrap();
    /// assert_eq!(x.checked_add(x), Q16::from_i32(4));
    /// assert_eq!(Q16::MAX.checked_add(Q16::ONE), None);
    /// ```
    pub fn checked_add(self, rhs: Self) -> Option<Self> {
        self.0.checked_add(rhs.0).map(Self)
    }

    /// Checked subtraction: exact, `None` on i64 raw overflow.
    ///
    /// Class: partially invertible; exact inverse of [`Fixed::checked_add`].
    pub fn checked_sub(self, rhs: Self) -> Option<Self> {
        self.0.checked_sub(rhs.0).map(Self)
    }

    /// Checked negation: exact, `None` only for `Fixed::MIN`.
    pub fn checked_neg(self) -> Option<Self> {
        self.0.checked_neg().map(Self)
    }

    /// Checked multiplication with floor rounding.
    ///
    /// Computes `r = floor(a_raw * b_raw / 2^FRAC)` in i128 and returns
    /// `None` if `r` does not fit in `i64`.
    ///
    /// Error bound (worst case, exact): the exact real product is
    /// `a_raw*b_raw / 2^(2*FRAC)` and the returned value is
    /// `r / 2^FRAC`, so `-2^-FRAC < err <= 0` (at most 1 ulp, rounded
    /// toward -infinity; 0 when the raw product is divisible by `2^FRAC`).
    /// Verified by `tests/fixed_bounds.rs::mul_floor_bound`.
    ///
    /// Class: lossy (many raw pairs map to one result); the domain guard
    /// (`None`) only covers representability, not the rounding.
    ///
    /// ```
    /// use phase_math::Q16;
    /// let x = Q16::from_i32(3).unwrap();
    /// assert_eq!(x.checked_mul(Q16::from_i32(7).unwrap()), Q16::from_i32(21));
    /// // 0.5 * 0.5 = 0.25 exactly (representable):
    /// let half = Q16::from_raw(1 << 15);
    /// assert_eq!(half.checked_mul(half), Some(Q16::from_raw(1 << 14)));
    /// // Floor rounding: 3 * 2^-48 collapses to 0 (error < 1 ulp):
    /// assert_eq!(Q16::from_raw(3).checked_mul(Q16::from_raw(1)), Some(Q16::ZERO));
    /// // Out of domain:
    /// assert_eq!(x.checked_mul(Q16::MAX), None);
    /// ```
    pub fn checked_mul(self, rhs: Self) -> Option<Self> {
        let p = (self.0 as i128) * (rhs.0 as i128);
        let r = p >> FRAC; // arithmetic shift = floor division
        if r < i64::MIN as i128 || r > i64::MAX as i128 {
            return None;
        }
        Some(Self(r as i64))
    }

    /// Saturating multiplication with floor rounding and a saturation flag.
    ///
    /// Same rounding as [`Fixed::checked_mul`]; an out-of-range result is
    /// clamped to [`Fixed::MAX`] / [`Fixed::MIN`] and the flag is set.
    /// Budget: flag `false` => exact within the multiplication bound above;
    /// flag `true` => no numeric bound, value is the clamp.
    pub fn saturating_mul(self, rhs: Self) -> Saturated<Self> {
        let p = (self.0 as i128) * (rhs.0 as i128);
        let r = p >> FRAC;
        let (v, saturated) = if r > i64::MAX as i128 {
            (i64::MAX, true)
        } else if r < i64::MIN as i128 {
            (i64::MIN, true)
        } else {
            (r as i64, false)
        };
        Saturated {
            value: Self(v),
            saturated,
        }
    }

    /// Checked division with floor rounding.
    ///
    /// Computes `r = floor(a_raw * 2^FRAC / b_raw)` in i128 (exact numerator)
    /// and returns `None` if `b == 0` or `r` does not fit in `i64`.
    ///
    /// Error bound (worst case, exact): the exact real quotient is
    /// `a_raw / b_raw` and the returned value is `r / 2^FRAC`, so
    /// `-2^-FRAC < err <= 0` (at most 1 ulp, toward -infinity; 0 when the
    /// division is exact). Verified by
    /// `tests/fixed_bounds.rs::div_floor_bound`.
    ///
    /// Note: this floors, while the legacy `q_div` truncates toward zero;
    /// both conventions are documented in `docs/math-spec.md`.
    ///
    /// Class: lossy; division by zero is a domain error (`None`).
    ///
    /// ```
    /// use phase_math::Q16;
    /// let a = Q16::from_i32(7).unwrap();
    /// assert_eq!(a.checked_div(Q16::from_i32(2).unwrap()), Some(Q16::from_raw(7 << 15))); // 3.5 exact
    /// assert_eq!(a.checked_div(Q16::ZERO), None); // div-by-0 is a domain error
    /// ```
    pub fn checked_div(self, rhs: Self) -> Option<Self> {
        if rhs.0 == 0 {
            return None;
        }
        debug_assert!(FRAC <= 63);
        let num = (self.0 as i128) << FRAC;
        let r = floor_div_i128(num, rhs.0 as i128);
        if r < i64::MIN as i128 || r > i64::MAX as i128 {
            return None;
        }
        Some(Self(r as i64))
    }

    /// Saturating addition with a saturation flag.
    ///
    /// Exact when the flag is `false`; clamped to `MAX`/`MIN` (flag `true`)
    /// otherwise. No numeric bound is offered when saturated.
    ///
    /// ```
    /// use phase_math::Q16;
    /// let s = Q16::MAX.saturating_add(Q16::ONE);
    /// assert!(s.saturated && s.value == Q16::MAX);
    /// let t = Q16::ONE.saturating_add(Q16::ONE);
    /// assert!(!t.saturated && t.value == Q16::from_raw(1 << 17));
    /// ```
    pub fn saturating_add(self, rhs: Self) -> Saturated<Self> {
        let value = Self(self.0.saturating_add(rhs.0));
        let saturated = value.0 != self.0.wrapping_add(rhs.0);
        Saturated { value, saturated }
    }

    /// Saturating subtraction with a saturation flag.
    ///
    /// Exact when the flag is `false`; clamped otherwise. See
    /// [`Fixed::saturating_add`] for the flag contract.
    pub fn saturating_sub(self, rhs: Self) -> Saturated<Self> {
        let value = Self(self.0.saturating_sub(rhs.0));
        let saturated = value.0 != self.0.wrapping_sub(rhs.0);
        Saturated { value, saturated }
    }

    /// Integer square root with floor rounding (requires `FRAC` even).
    ///
    /// Computes `r = isqrt(a_raw) * 2^(FRAC/2)`; returns `None` if the value
    /// is negative or `FRAC` is odd.
    ///
    /// Error bound (worst case, exact): the exact real result is
    /// `sqrt(a_raw) / 2^(FRAC/2)` and the returned value is
    /// `isqrt(a_raw) / 2^(FRAC/2)`, so `-2^-(FRAC/2) < err <= 0`
    /// (for `Q16`: `-1/256 < err <= 0`, i.e. at most 256 ulp — squaring
    /// halves the precision). Verified by `tests/fixed_bounds.rs::sqrt_bound`.
    ///
    /// Class: lossy (strictly monotone but many-to-one).
    ///
    /// ```
    /// use phase_math::Q16;
    /// assert_eq!(Q16::from_i32(16).unwrap().sqrt(), Q16::from_i32(4));
    /// // floor(sqrt(2)) * 2^-8: -2^-8 < err <= 0
    /// assert_eq!(Q16::from_raw(2).sqrt(), Some(Q16::from_raw(256)));
    /// ```
    pub fn sqrt(self) -> Option<Self> {
        if self.0 < 0 || FRAC % 2 != 0 {
            return None;
        }
        Some(Self((isqrt_u64(self.0 as u64) << (FRAC / 2)) as i64))
    }

    /// `2^e` in fixed-point.
    ///
    /// Exact (error 0) for `-FRAC <= e <= 62 - FRAC` (raw `2^(FRAC+e)` must
    /// fit `i64`). For `e <= -FRAC - 1` the exact raw `2^(FRAC+e)` is below
    /// 1 ulp and the result is `0` with error `2^e < 2^-FRAC`. Returns
    /// `None` for `e > 62 - FRAC` (unrepresentable).
    ///
    /// Class: partially invertible — [`Fixed::log2`] is the exact inverse on
    /// `e ∈ [-FRAC, 62 - FRAC]`; all `e <= -FRAC - 1` collapse to `0`.
    ///
    /// ```
    /// use phase_math::Q16;
    /// assert_eq!(Q16::pow2(2), Q16::from_i32(4));       // 4.0 exact
    /// assert_eq!(Q16::pow2(-16), Some(Q16::from_raw(1))); // 2^-16 = 1 ulp, exact
    /// assert_eq!(Q16::pow2(-17), Some(Q16::ZERO));      // below resolution: 0, err = 2^-17
    /// assert_eq!(Q16::pow2(47), None);                  // 2^63 does not fit
    /// ```
    pub fn pow2(e: i32) -> Option<Self> {
        let sh = FRAC as i64 + e as i64; // raw = 2^sh
        if sh < 0 {
            return Some(Self(0));
        }
        if sh > 62 {
            return None; // 2^63 does not fit in i64
        }
        Some(Self(1i64 << sh))
    }

    /// `floor(log2(value))` as `i32`.
    ///
    /// Returns `None` for values `<= 0` (domain). The result is the exact
    /// floor of the real `log2`: error `-1 < err_log2 <= 0` in bits. This is
    /// a many-to-one map (lossy class); it is the exact inverse of
    /// [`Fixed::pow2`] on `e ∈ [-FRAC, 62 - FRAC]`.
    ///
    /// ```
    /// use phase_math::Q16;
    /// assert_eq!(Q16::from_i32(8).unwrap().log2(), Some(3));
    /// assert_eq!(Q16::from_raw(3).log2(), Some(-15)); // floor(log2(3) - 16)
    /// assert_eq!(Q16::ZERO.log2(), None);
    /// ```
    pub fn log2(self) -> Option<i32> {
        if self.0 <= 0 {
            return None;
        }
        Some(63 - self.0.leading_zeros() as i32 - FRAC as i32)
    }
}

impl<const FRAC: u32> fmt::Display for Fixed<FRAC> {
    /// Exact decimal rendering: the value `raw / 2^FRAC` has a finite
    /// decimal expansion with at most `FRAC` fractional digits; trailing
    /// zeros are stripped. No rounding happens in `Display`.
    ///
    /// ```
    /// use phase_math::Q16;
    /// assert_eq!(format!("{}", Q16::from_raw(32768)), "0.5");
    /// assert_eq!(format!("{}", Q16::from_raw(-196609)), "-3.0000152587890625");
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        debug_assert!(FRAC <= 63);
        let neg = self.0 < 0;
        let n = self.0.unsigned_abs() as u128;
        let mask = (1u128 << FRAC) - 1;
        let int_part = n >> FRAC;
        let mut frac = n & mask;
        let mut digits = [0u8; 64]; // FRAC <= 63
        let mut len = 0usize;
        for _ in 0..FRAC {
            frac *= 10;
            digits[len] = b'0' + (frac >> FRAC) as u8;
            len += 1;
            frac &= mask;
        }
        while len > 0 && digits[len - 1] == b'0' {
            len -= 1;
        }
        if neg && (int_part != 0 || len > 0) {
            write!(f, "-")?;
        }
        write!(f, "{}", int_part)?;
        if len > 0 {
            write!(f, ".")?;
            for d in &digits[..len] {
                write!(f, "{}", *d as char)?;
            }
        }
        Ok(())
    }
}

/// Floor division for i128 (`b != 0`): `q = floor(a / b)`.
///
/// Property (verified by `tests/fixed_bounds.rs::div_floor_bound` and the
/// interval tests): the remainder `a - q*b` has the sign of `b` (or is
/// zero) and `|a - q*b| < |b|`. For the represented real values this means
/// the computed quotient never exceeds the exact one (one-sided error,
/// toward -infinity), regardless of the divisor's sign.
pub(crate) const fn floor_div_i128(a: i128, b: i128) -> i128 {
    let q = a / b;
    let r = a % b;
    if r != 0 && (a < 0) != (b < 0) {
        q - 1
    } else {
        q
    }
}

/// Ceiling division for i128 (`b != 0`): `q = ceil(a / b)`.
///
/// Property: the remainder `q*b - a` has the sign of `b` (or is zero) and
/// `|q*b - a| < |b|` — the computed quotient never falls below the exact one.
pub(crate) const fn ceil_div_i128(a: i128, b: i128) -> i128 {
    let q = a / b;
    let r = a % b;
    if r != 0 && (a < 0) == (b < 0) {
        q + 1
    } else {
        q
    }
}

/// Integer square root (floor), Newton's method, 64-bit.
///
/// Overflow-free for the whole `u64` range: the initial guess is `n/2 + 1`
/// (never `n + 1`) and every Newton step stays below `2^64`.
pub(crate) fn isqrt_u64(n: u64) -> u64 {
    if n < 2 {
        return n;
    }
    let mut x = n;
    let mut y = x / 2 + (x & 1); // == (x + 1) / 2 without the overflow
    while y < x {
        x = y;
        y = (x + n / x) / 2;
    }
    x
}

/// Integer square root (ceiling): smallest `s` with `s*s >= n`.
pub(crate) fn ceil_sqrt_u64(n: u64) -> u64 {
    let s = isqrt_u64(n);
    if s.saturating_mul(s) < n {
        s + 1
    } else {
        s
    }
}

/// `Fixed` with 16 fractional bits — the classic "Q16.16" layout used by
/// the legacy `q_*` functions (`type Q = i64` raw).
pub type Q16 = Fixed<16>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn floor_ceil_div_properties() {
        let mut cases: Vec<(i128, i128)> = vec![
            (5, 2),
            (-5, 2),
            (5, -2),
            (-5, -2),
            (7, 3),
            (-7, 3),
            (0, 5),
            (i64::MIN as i128, 1),
            (i64::MIN as i128, -1),
            (i64::MAX as i128, 2),
            (1 << 79, 1 << 16),
            (-(1 << 79), 1 << 16),
        ];
        // deterministic pseudo-random sweep
        let mut s = 0x1234_5678_9ABC_DEF0u64;
        let mut next = move || {
            s ^= s >> 12;
            s ^= s << 25;
            s ^= s >> 27;
            s.wrapping_mul(0x2545_F491_4F6C_DD1D)
        };
        for _ in 0..20_000 {
            let a = next() as i64 as i128 * (next() % 4 + 1) as i128;
            let mut b = next() as i64 as i128;
            if b == 0 {
                b = 1;
            }
            cases.push((a, b));
        }
        for &(a, b) in &cases {
            let fq = floor_div_i128(a, b);
            let rem = a - fq * b;
            // floor: remainder has the sign of b, |rem| < |b|
            assert!(rem == 0 || (rem > 0) == (b > 0), "floor {a}/{b} -> {fq}");
            assert!(rem.abs() < b.abs(), "floor {a}/{b} -> {fq}");
            let cq = ceil_div_i128(a, b);
            let rem2 = cq * b - a;
            // ceil: same sign property on q*b - a
            assert!(rem2 == 0 || (rem2 > 0) == (b > 0), "ceil {a}/{b} -> {cq}");
            assert!(rem2.abs() < b.abs(), "ceil {a}/{b} -> {cq}");
            assert!(fq <= cq);
        }
    }

    #[test]
    fn isqrt_and_ceil_sqrt() {
        for n in 0u64..10_000 {
            let s = isqrt_u64(n);
            assert!(
                (s as u128) * (s as u128) <= n as u128
                    && (n as u128) < ((s + 1) as u128) * ((s + 1) as u128),
                "isqrt({n}) = {s}"
            );
            let c = ceil_sqrt_u64(n);
            assert!(
                (c as u128) * (c as u128) >= n as u128
                    && (c == 0 || ((c - 1) as u128) * ((c - 1) as u128) < n as u128),
                "ceil_sqrt({n}) = {c}"
            );
        }
        // 32-bit boundary: (2^32 - 1)^2 < 2^64, isqrt must not overflow
        let big = u64::MAX;
        assert_eq!(isqrt_u64(big), u32::MAX as u64);
        assert_eq!(ceil_sqrt_u64(big), 1u64 << 32);
        for _ in 0..50_000 {
            let n = 0x1234_5678_9ABC_DEF0u64
                .wrapping_mul(0x2545_F491_4F6C_DD1D)
                .wrapping_add(1);
            let s = isqrt_u64(n);
            assert!(
                (s as u128) * (s as u128) <= n as u128
                    && (n as u128) < ((s + 1) as u128) * ((s + 1) as u128)
            );
        }
    }

    #[test]
    fn display_exact_decimal() {
        assert_eq!(format!("{}", Q16::from_raw(0)), "0");
        assert_eq!(format!("{}", Q16::from_raw(1 << 16)), "1");
        assert_eq!(format!("{}", Q16::from_raw(-(1 << 16))), "-1");
        assert_eq!(format!("{}", Q16::from_raw(1)), "0.0000152587890625");
        assert_eq!(format!("{}", Q16::from_raw(-1)), "-0.0000152587890625");
        assert_eq!(format!("{}", Q16::from_raw(65536 + 32768)), "1.5");
        assert_eq!(format!("{}", Fixed::<8>::from_raw(-256)), "-1");
        assert_eq!(format!("{}", Fixed::<8>::from_raw(-255)), "-0.99609375");
        assert_eq!(format!("{}", Fixed::<1>::from_raw(1)), "0.5");
    }

    #[test]
    fn q16_constants_and_ordering() {
        assert_eq!(Q16::ONE.raw(), 1 << 16);
        assert_eq!(Q16::ZERO.raw(), 0);
        assert!(Q16::MIN < Q16::ZERO && Q16::ZERO < Q16::ONE && Q16::ONE < Q16::MAX);
        assert_eq!(Q16::default(), Q16::ZERO);
    }

    #[test]
    fn sqrt_domain_and_odd_frac() {
        assert_eq!(Q16::from_raw(-1).sqrt(), None);
        assert_eq!(Fixed::<15>::from_raw(100).sqrt(), None); // odd FRAC -> None
        assert_eq!(Fixed::<0>::from_raw(9).sqrt(), Some(Fixed::<0>::from_raw(3)));
    }

    #[test]
    fn from_i32_domain() {
        assert_eq!(Fixed::<32>::from_i32(i32::MAX), Some(Fixed::<32>::from_raw((i32::MAX as i64) << 32)));
        assert_eq!(Fixed::<33>::from_i32(i32::MIN), None); // -2^63 << 1 overflows
    }
}
