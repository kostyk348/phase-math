//! Interval arithmetic over [`Fixed`] with outward rounding.
//!
//! An [`Interval`]`<FRAC>` is a closed pair `lo <= hi` of fixed-point
//! endpoints. Every operation is **outward-rounded**: for any real inputs
//! within the operand intervals, the exact real result is guaranteed to lie
//! inside the returned interval. Soundness is never traded for
//! representability — when an endpoint cannot be represented in `i64` raw,
//! the operation returns `None` instead of inflating toward infinity (this
//! crate has no infinities; see the honesty notes in `docs/math-spec.md`).
//!
//! Tightness (also verified by tests): for multiplication and division the
//! endpoints are exactly the true extreme corners rounded outward (floor
//! for `lo`, ceiling for `hi`), so each endpoint is within 1 ulp
//! (`2^-FRAC`) of the exact extreme, and the interval width exceeds the
//! exact width by at most 2 ulp. Addition and subtraction are exact and
//! tight (zero inflation).

use crate::fixed::{ceil_div_i128, ceil_sqrt_u64, floor_div_i128, isqrt_u64, Fixed};

/// A closed interval `[lo, hi]` over `Fixed<FRAC>`, maintained with
/// outward rounding: every real result of an operation on values inside
/// the operands is contained in the returned interval.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Interval<const FRAC: u32> {
    /// Lower endpoint (inclusive).
    pub lo: Fixed<FRAC>,
    /// Upper endpoint (inclusive).
    pub hi: Fixed<FRAC>,
}

impl<const FRAC: u32> Interval<FRAC> {
    /// Builds an interval. Returns `None` if `lo > hi` (empty intervals are
    /// not representable in this type).
    ///
    /// ```
    /// use phase_math::{Interval, Q16};
    /// assert!(Interval::new(Q16::ZERO, Q16::ONE).is_some());
    /// assert!(Interval::new(Q16::ONE, Q16::ZERO).is_none());
    /// ```
    pub fn new(lo: Fixed<FRAC>, hi: Fixed<FRAC>) -> Option<Self> {
        if lo <= hi {
            Some(Self { lo, hi })
        } else {
            None
        }
    }

    /// A degenerate interval containing exactly one value.
    pub fn point(x: Fixed<FRAC>) -> Self {
        Self { lo: x, hi: x }
    }

    /// `true` iff `x` lies inside `[lo, hi]` (endpoints inclusive).
    pub fn contains(self, x: Fixed<FRAC>) -> bool {
        self.lo <= x && x <= self.hi
    }

    /// `true` iff the interval contains the exact value 0. Division by such
    /// an interval is a domain error (`None`), because the exact quotient
    /// may be unbounded.
    pub fn contains_zero(self) -> bool {
        self.lo.raw() <= 0 && 0 <= self.hi.raw()
    }

    /// Width `hi - lo`; exact, `None` on raw overflow.
    pub fn width(self) -> Option<Fixed<FRAC>> {
        self.hi.checked_sub(self.lo)
    }

    /// Interval addition: exact and tight (zero inflation).
    ///
    /// `lo = a.lo + b.lo`, `hi = a.hi + b.hi`, both computed exactly;
    /// returns `None` if an endpoint overflows `i64` raw.
    ///
    /// Guarantee: for any real `x ∈ [a.lo, a.hi]`, `y ∈ [b.lo, b.hi]`, the
    /// exact sum `x + y` lies in the returned interval.
    ///
    /// ```
    /// use phase_math::{Interval, Q16};
    /// let a = Interval::new(Q16::from_raw(0), Q16::ONE).unwrap();
    /// let s = a.add(a).unwrap(); // [0, 2] exactly
    /// assert_eq!(s.lo, Q16::ZERO);
    /// assert_eq!(s.hi, Q16::from_raw(2 << 16));
    /// ```
    pub fn add(self, rhs: Self) -> Option<Self> {
        Some(Self {
            lo: self.lo.checked_add(rhs.lo)?,
            hi: self.hi.checked_add(rhs.hi)?,
        })
    }

    /// Interval subtraction: exact and tight (zero inflation).
    ///
    /// `lo = a.lo - b.hi`, `hi = a.hi - b.lo`; `None` on raw overflow.
    pub fn sub(self, rhs: Self) -> Option<Self> {
        Some(Self {
            lo: self.lo.checked_sub(rhs.hi)?,
            hi: self.hi.checked_sub(rhs.lo)?,
        })
    }

    /// Interval multiplication with outward rounding.
    ///
    /// All four corner products `a_i * b_j` are computed exactly in i128;
    /// `lo` is the minimum of their floor-roundings and `hi` the maximum of
    /// their ceiling-roundings. Returns `None` if an endpoint is not
    /// representable.
    ///
    /// Guarantee: the interval contains every exact real product
    /// `x * y` with `x` in `a`, `y` in `b` (products are bilinear, so the
    /// extremes are attained at corners). Tightness: each endpoint is the
    /// outward rounding of the exact extreme corner, hence within 1 ulp of
    /// it; width exceeds the exact width by at most 2 ulp. Verified by
    /// `tests/interval_outward.rs::mul_containment_and_exact_outward_ends`.
    ///
    /// Class: per-endpoint lossy rounding (outward), containment exact.
    ///
    /// ```
    /// use phase_math::{Interval, Q16};
    /// let a = Interval::new(Q16::from_raw(3 * (1 << 16)), Q16::from_raw(4 * (1 << 16))).unwrap();
    /// let b = Interval::new(Q16::from_raw(-1), Q16::from_raw(2)).unwrap(); // [~0, 2 ulp]
    /// let m = a.mul(b).unwrap();
    /// assert!(m.lo.raw() <= 0 && 0 <= m.hi.raw()); // exact products in [-0, 8 ulp]
    /// ```
    pub fn mul(self, rhs: Self) -> Option<Self> {
        let (al, ah) = (self.lo.raw() as i128, self.hi.raw() as i128);
        let (bl, bh) = (rhs.lo.raw() as i128, rhs.hi.raw() as i128);
        let mut lo = i128::MAX;
        let mut hi = i128::MIN;
        for &a in &[al, ah] {
            for &b in &[bl, bh] {
                let p = a * b;
                lo = lo.min(p >> FRAC); // floor(p / 2^FRAC)
                hi = hi.max(-((-p) >> FRAC)); // ceil(p / 2^FRAC)
            }
        }
        if lo < i64::MIN as i128 || hi > i64::MAX as i128 {
            return None;
        }
        Some(Self {
            lo: Fixed::from_raw(lo as i64),
            hi: Fixed::from_raw(hi as i64),
        })
    }

    /// Interval division with outward rounding.
    ///
    /// Returns `None` if the divisor interval contains 0 (the exact
    /// quotient may be unbounded) or an endpoint is not representable.
    /// Otherwise `lo` is the minimum of the exact corner quotients rounded
    /// down and `hi` the maximum rounded up, each quotient computed with an
    /// exact i128 numerator `a_i * 2^FRAC`.
    ///
    /// Guarantee: contains every exact real quotient `x / y` with `x` in
    /// `a`, `y` in `b`; endpoints are the outward roundings of the exact
    /// extreme corners (within 1 ulp; width +2 ulp at most). Verified by
    /// `tests/interval_outward.rs::div_containment_and_exact_outward_ends`.
    ///
    /// ```
    /// use phase_math::{Interval, Q16};
    /// let a = Interval::point(Q16::from_i32(7).unwrap());
    /// let b = Interval::new(Q16::from_raw(1 << 15), Q16::from_raw(1 << 16)).unwrap(); // [0.5, 1]
    /// let q = a.div(b).unwrap(); // 7 / [0.5, 1] = [7, 14]
    /// assert_eq!(q.lo, Q16::from_raw(7 << 16));
    /// assert_eq!(q.hi, Q16::from_raw(14 << 16));
    /// ```
    pub fn div(self, rhs: Self) -> Option<Self> {
        if rhs.contains_zero() {
            return None;
        }
        debug_assert!(FRAC <= 63);
        let (al, ah) = (self.lo.raw() as i128, self.hi.raw() as i128);
        let (bl, bh) = (rhs.lo.raw() as i128, rhs.hi.raw() as i128);
        let mut lo = i128::MAX;
        let mut hi = i128::MIN;
        for &a in &[al, ah] {
            for &b in &[bl, bh] {
                let num = a << FRAC;
                lo = lo.min(floor_div_i128(num, b));
                hi = hi.max(ceil_div_i128(num, b));
            }
        }
        if lo < i64::MIN as i128 || hi > i64::MAX as i128 {
            return None;
        }
        Some(Self {
            lo: Fixed::from_raw(lo as i64),
            hi: Fixed::from_raw(hi as i64),
        })
    }

    /// Interval square root with outward rounding (requires `FRAC` even and
    /// `lo >= 0`).
    ///
    /// `lo = floor_sqrt(a.lo_raw) * 2^(FRAC/2)`, `hi = ceil_sqrt(a.hi_raw) *
    /// 2^(FRAC/2)`; `None` if the interval reaches into negative values or
    /// `FRAC` is odd (sqrt is not real there).
    ///
    /// Guarantee: contains every exact real `sqrt(x)`, `x` in the operand.
    /// Each endpoint is within `2^-(FRAC/2)` of the exact extreme (the same
    /// bound as scalar [`Fixed::sqrt`]). Verified by
    /// `tests/interval_outward.rs::sqrt_containment`.
    pub fn sqrt(self) -> Option<Self> {
        if self.lo.raw() < 0 || FRAC % 2 != 0 {
            return None;
        }
        let lo = (isqrt_u64(self.lo.raw() as u64) << (FRAC / 2)) as i64;
        let hi = (ceil_sqrt_u64(self.hi.raw() as u64) << (FRAC / 2)) as i64;
        Some(Self {
            lo: Fixed::from_raw(lo),
            hi: Fixed::from_raw(hi),
        })
    }
}
