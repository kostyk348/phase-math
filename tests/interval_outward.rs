//! Interval containment and tightness tests against exact i128 rational
//! references (docs/math-spec.md, Family E).
//!
//! For every operation: (1) all exact real corner results must lie inside
//! the returned interval (containment), and (2) the endpoints must be
//! exactly the outward roundings of the true extreme corners (tightness —
//! hence within 1 ulp of the exact extreme, width +2 ulp at most).

use phase_math::{Fixed, Interval, Q16};

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
}

fn ref_floor_div(a: i128, b: i128) -> i128 {
    let q = a / b;
    let r = a % b;
    if r != 0 && (a < 0) != (b < 0) {
        q - 1
    } else {
        q
    }
}

fn ref_ceil_div(a: i128, b: i128) -> i128 {
    let q = a / b;
    let r = a % b;
    if r != 0 && (a < 0) == (b < 0) {
        q + 1
    } else {
        q
    }
}

/// Independent reference: floor sqrt via binary search (no Newton, no f64).
fn isqrt_ref(n: u64) -> u64 {
    let mut lo: u64 = 0;
    let mut hi: u64 = 1 << 32;
    while lo + 1 < hi {
        let mid = (lo + hi) / 2;
        if mid * mid <= n {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo
}

/// Independent reference: ceiling sqrt via floor sqrt + correction.
fn ceil_sqrt_ref(n: u64) -> u64 {
    let lo = isqrt_ref(n);
    if lo * lo < n {
        lo + 1
    } else {
        lo
    }
}

fn rand_interval(rng: &mut Rng, span: i64) -> Interval<16> {
    let half = 2u64 * (span as u64) + 1;
    let a = (rng.next() % half) as i64 - span;
    let b = (rng.next() % half) as i64 - span;
    Interval::new(Q16::from_raw(a.min(b)), Q16::from_raw(a.max(b))).unwrap()
}

fn rand_full_interval(rng: &mut Rng) -> Interval<16> {
    let a = rng.next() as i64;
    let b = rng.next() as i64;
    Interval::new(Q16::from_raw(a.min(b)), Q16::from_raw(a.max(b))).unwrap()
}

#[test]
fn mul_containment_and_exact_outward_ends() {
    let mut rng = Rng(0x17E_8A11);
    for _ in 0..50_000 {
        let a = rand_interval(&mut rng, 1 << 30);
        let b = rand_interval(&mut rng, 1 << 30);
        let m = a.mul(b).expect("corners of +-2^30 fit for F=16");
        let corners: [i128; 4] = [
            (a.lo.raw() as i128) * (b.lo.raw() as i128),
            (a.lo.raw() as i128) * (b.hi.raw() as i128),
            (a.hi.raw() as i128) * (b.lo.raw() as i128),
            (a.hi.raw() as i128) * (b.hi.raw() as i128),
        ];
        let p_min = *corners.iter().min().unwrap();
        let p_max = *corners.iter().max().unwrap();
        // 1) containment: every exact corner result inside [lo, hi]
        //    lo/2^16 <= p/2^32 <= hi/2^16  <=>  lo*2^16 <= p <= hi*2^16
        for p in corners {
            assert!((m.lo.raw() as i128) * (1 << 16) <= p, "lo containment");
            assert!(p <= (m.hi.raw() as i128) * (1 << 16), "hi containment");
        }
        // 2) exact outward ends: endpoints are the outward roundings of the
        //    true extreme corners (floor for lo, ceil for hi)
        assert_eq!(m.lo.raw() as i128, p_min >> 16);
        assert_eq!(m.hi.raw() as i128, ref_ceil_div(p_max, 1 << 16));
        // 3) width <= exact width + 2 ulp
        let w = (m.hi.raw() - m.lo.raw()) as i128;
        assert!(w * (1 << 16) <= p_max - p_min + 2 * (1 << 16));
    }
}

#[test]
fn div_containment_and_exact_outward_ends() {
    let mut rng = Rng(0x17E_8A22);
    for case in 0..50_000u64 {
        let a = rand_interval(&mut rng, 1 << 30);
        // strictly positive / strictly negative divisor interval
        let (u, v) = (
            1 + (rng.next() % (1 << 29)) as i64,
            1 + (rng.next() % (1 << 29)) as i64,
        );
        let b = if case % 2 == 0 {
            Interval::new(Q16::from_raw(u.min(v)), Q16::from_raw(u.max(v))).unwrap()
        } else {
            Interval::new(Q16::from_raw(-u.max(v)), Q16::from_raw(-u.min(v))).unwrap()
        };
        let q = a.div(b).expect("quotients fit for these spans");
        // exact corner quotients: num = ai * 2^16 over den = bj (i128 exact)
        let mut corners: Vec<(i128, i128)> = Vec::with_capacity(4);
        for &ai in &[a.lo.raw() as i128, a.hi.raw() as i128] {
            for &bj in &[b.lo.raw() as i128, b.hi.raw() as i128] {
                corners.push((ai << 16, bj));
            }
        }
        // 1) containment with sign-aware cross multiplication:
        //    lo <= num/den  <=>  (den > 0) ? lo*den <= num : num <= lo*den
        for &(num, den) in &corners {
            let lo = q.lo.raw() as i128;
            let hi = q.hi.raw() as i128;
            let inside_lo = if den > 0 { lo * den <= num } else { num <= lo * den };
            let inside_hi = if den > 0 { num <= hi * den } else { hi * den <= num };
            assert!(inside_lo, "lo containment {num}/{den}");
            assert!(inside_hi, "hi containment {num}/{den}");
        }
        // 2) exact outward ends: lo is the floor of the true minimum corner,
        //    hi the ceil of the true maximum corner (compare fractions by
        //    cross multiplication; all dens share a sign)
        let min_c = *corners
            .iter()
            .reduce(|c1, c2| {
                if c1.0 * c2.1 < c2.0 * c1.1 {
                    c1
                } else {
                    c2
                }
            })
            .unwrap();
        let max_c = *corners
            .iter()
            .reduce(|c1, c2| {
                if c1.0 * c2.1 > c2.0 * c1.1 {
                    c1
                } else {
                    c2
                }
            })
            .unwrap();
        assert_eq!(q.lo.raw() as i128, ref_floor_div(min_c.0, min_c.1));
        assert_eq!(q.hi.raw() as i128, ref_ceil_div(max_c.0, max_c.1));
    }
}

#[test]
fn add_sub_exact_and_tight() {
    let mut rng = Rng(0x17E_8A33);
    for _ in 0..25_000 {
        let a = rand_interval(&mut rng, 1 << 40);
        let b = rand_interval(&mut rng, 1 << 40);
        // add: endpoints are the exact sums; None iff an endpoint overflows
        match a.add(b) {
            Some(s) => {
                assert_eq!(s.lo.raw(), a.lo.raw() + b.lo.raw());
                assert_eq!(s.hi.raw(), a.hi.raw() + b.hi.raw());
                // zero inflation: width == exact width
                assert_eq!(s.width().unwrap().raw(), (a.hi.raw() - a.lo.raw()) + (b.hi.raw() - b.lo.raw()));
            }
            None => assert!(
                a.lo.raw().checked_add(b.lo.raw()).is_none()
                    || a.hi.raw().checked_add(b.hi.raw()).is_none()
            ),
        }
        // sub: lo = a.lo - b.hi, hi = a.hi - b.lo (exact)
        match a.sub(b) {
            Some(s) => {
                assert_eq!(s.lo.raw(), a.lo.raw() - b.hi.raw());
                assert_eq!(s.hi.raw(), a.hi.raw() - b.lo.raw());
            }
            None => assert!(
                a.lo.raw().checked_sub(b.hi.raw()).is_none()
                    || a.hi.raw().checked_sub(b.lo.raw()).is_none()
            ),
        }
    }
    // full-range endpoints exercise the overflow -> None paths
    for _ in 0..25_000 {
        let a = rand_full_interval(&mut rng);
        let b = rand_full_interval(&mut rng);
        match a.add(b) {
            Some(s) => {
                assert_eq!(s.lo.raw(), a.lo.raw().wrapping_add(b.lo.raw()));
                assert_eq!(s.hi.raw(), a.hi.raw().wrapping_add(b.hi.raw()));
            }
            None => assert!(
                a.lo.raw().checked_add(b.lo.raw()).is_none()
                    || a.hi.raw().checked_add(b.hi.raw()).is_none()
            ),
        }
    }
}

#[test]
fn sqrt_containment_and_outward_ends() {
    let mut rng = Rng(0x17E_8A44);
    for _ in 0..50_000 {
        let lo = (rng.next() >> 4) as i64; // non-negative raws
        let hi = lo + (rng.next() % (1 << 24)) as i64;
        let a = Interval::new(Q16::from_raw(lo), Q16::from_raw(hi)).unwrap();
        let s = a.sqrt().unwrap();
        // exact outward ends by construction:
        let l = isqrt_ref(lo as u64);
        let h = ceil_sqrt_ref(hi as u64);
        assert_eq!(s.lo.raw() as u64, l << 8);
        assert_eq!(s.hi.raw() as u64, h << 8);
        // containment: floor(lo)^2 <= lo and ceil(hi)^2 >= hi
        assert!(l * l <= (lo as u64) && (lo as u64) < (l + 1) * (l + 1));
        assert!(h * h >= (hi as u64) && (h == 0 || (h - 1) * (h - 1) < (hi as u64)));
    }
    // domain: negative lo -> None
    assert_eq!(
        Interval::new(Q16::from_raw(-1), Q16::ONE).unwrap().sqrt(),
        None
    );
    // odd FRAC -> None
    assert!(
        Interval::new(Fixed::<15>::from_raw(0), Fixed::<15>::from_raw(9))
            .unwrap()
            .sqrt()
            .is_none()
    );
    // the whole representable range stays representable after sqrt
    let full = Interval::new(Q16::ZERO, Q16::MAX).unwrap().sqrt().unwrap();
    assert_eq!(full.lo, Q16::ZERO);
    assert!(full.hi.raw() > 0);
    // sqrt of [0, 2^62] contains 2^31 exactly
    let s = Interval::new(Q16::ZERO, Q16::from_raw(1i64 << 62))
        .unwrap()
        .sqrt()
        .unwrap();
    assert!(s.contains(Q16::from_raw(1i64 << 31)));
}

#[test]
fn domain_errors_and_basics() {
    let zero_one = Interval::new(Q16::from_raw(-1), Q16::ONE).unwrap();
    // dividing BY an interval containing 0 -> None (quotient may be unbounded)
    assert_eq!(Interval::point(Q16::ONE).div(zero_one), None);
    // but an interval containing 0 in the NUMERATOR is fine: [-1,1]/1 = [-1,1]
    assert_eq!(zero_one.div(Interval::point(Q16::ONE)), Some(zero_one));
    assert!(zero_one.contains_zero());
    assert!(!Interval::new(Q16::from_raw(1), Q16::from_raw(2))
        .unwrap()
        .contains_zero());
    // arithmetic overflow -> None (no infinities: refuse instead of inflating)
    assert_eq!(Interval::point(Q16::MAX).add(Interval::point(Q16::ONE)), None);
    assert_eq!(Interval::point(Q16::MIN).sub(Interval::point(Q16::ONE)), None);
    assert_eq!(Interval::new(Q16::MIN, Q16::MAX).unwrap().mul(Interval::new(Q16::MIN, Q16::MAX).unwrap()), None);
    // exact quotient unrepresentable: MAX_raw / 1_raw = (2^63-1) needs raw 2^79
    assert_eq!(
        Interval::point(Q16::MAX).div(Interval::point(Q16::from_raw(1))),
        None
    );
    // constructor: inverted bounds rejected
    assert_eq!(Interval::new(Q16::ONE, Q16::ZERO), None);
    // contains / point / width
    let i = Interval::new(Q16::from_raw(-65536), Q16::from_raw(131072)).unwrap();
    assert!(i.contains(Q16::ZERO) && i.contains(Q16::from_raw(131072)));
    assert!(!i.contains(Q16::from_raw(131073)));
    assert_eq!(Interval::point(Q16::ONE).width(), Some(Q16::ZERO));
    assert_eq!(i.width(), Some(Q16::from_raw(196608)));
}

#[test]
fn interval_display_and_equality_smoke() {
    // endpoints carry Fixed's exact Display
    let i = Interval::new(Q16::from_raw(32768), Q16::from_raw(65536)).unwrap();
    assert_eq!(format!("{}..{}", i.lo, i.hi), "0.5..1");
    // derived traits work
    let j = Interval {
        lo: Q16::from_raw(32768),
        hi: Q16::from_raw(65536),
    };
    assert_eq!(i, j);
    assert_ne!(i, Interval::point(Q16::from_raw(32768)));
}
