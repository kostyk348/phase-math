//! Error-bound tests against exact i128 references.
//!
//! Every worst-case bound claimed in docs/math-spec.md for the `Fixed<FRAC>`
//! families is checked here on edge cases plus large random samples. The
//! references use exact i128 rational arithmetic — no floats anywhere.

use phase_math::{Fixed, Q16};

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    fn i64(&mut self) -> i64 {
        self.next() as i64
    }
}

const EDGES: [i64; 10] = [
    0,
    1,
    -1,
    2,
    -2,
    i64::MAX,
    i64::MIN,
    i64::MAX - 1,
    i64::MIN + 1,
    1 << 16,
];

/// Independent reference: floor division (binary-search-free, direct).
fn floor_div_ref(a: i128, b: i128) -> i128 {
    let q = a / b;
    let r = a % b;
    if r != 0 && (a < 0) != (b < 0) {
        q - 1
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

fn check_mul_f16(a: i64, b: i64) {
    let p = (a as i128) * (b as i128);
    let expected = p >> 16; // exact floor reference in i128
    match Q16::from_raw(a).checked_mul(Q16::from_raw(b)) {
        Some(c) => {
            assert_eq!(c.raw() as i128, expected, "mul {a} * {b}");
            // Floor property: 0 <= p - c*2^16 < 2^16, hence the worst-case
            // bound of docs/math-spec.md: -2^-16 < err <= 0 (<= 1 ulp, down).
            let num_err = p - ((c.raw() as i128) << 16);
            assert!(
                (0..(1i128 << 16)).contains(&num_err),
                "mul {a} * {b}: numerator error {num_err} outside [0, 2^16)"
            );
        }
        None => assert!(
            expected < i64::MIN as i128 || expected > i64::MAX as i128,
            "None must mean out-of-domain: {a} * {b}"
        ),
    }
}

fn check_div_f16(a: i64, b: i64) {
    if b == 0 {
        assert_eq!(Q16::from_raw(a).checked_div(Q16::ZERO), None);
        return;
    }
    let num = (a as i128) << 16;
    let expected = floor_div_ref(num, b as i128);
    match Q16::from_raw(a).checked_div(Q16::from_raw(b)) {
        Some(c) => {
            assert_eq!(c.raw() as i128, expected, "div {a} / {b}");
            // Floor property: rem = num - c*b has the sign of b (or is 0)
            // and |rem| < |b|; dividing by a negative b flips the remainder
            // back, so the worst-case bound -2^-16 < err <= 0 holds for any
            // b != 0 (docs/math-spec.md, Family C).
            let rem = num - (c.raw() as i128) * (b as i128);
            assert!(
                rem == 0 || (rem > 0) == (b > 0),
                "div {a} / {b}: remainder {rem} has wrong sign"
            );
            assert!(
                rem.abs() < (b as i128).abs(),
                "div {a} / {b}: |{rem}| >= |b|"
            );
        }
        None => assert!(
            expected < i64::MIN as i128 || expected > i64::MAX as i128,
            "None must mean out-of-domain: {a} / {b}"
        ),
    }
}

#[test]
fn mul_floor_bound_f16() {
    for a in EDGES {
        for b in EDGES {
            check_mul_f16(a, b);
        }
    }
    let mut rng = Rng(0xB00D_0001);
    let mut max_num_err: i128 = 0;
    for _ in 0..200_000 {
        let (a, b) = (rng.i64(), rng.i64());
        let p = (a as i128) * (b as i128);
        let expected = p >> 16;
        if let Some(c) = Q16::from_raw(a).checked_mul(Q16::from_raw(b)) {
            assert_eq!(c.raw() as i128, expected);
            let num_err = p - ((c.raw() as i128) << 16);
            assert!((0..(1i128 << 16)).contains(&num_err));
            max_num_err = max_num_err.max(num_err);
        } else {
            assert!(expected < i64::MIN as i128 || expected > i64::MAX as i128);
        }
    }
    // The measured maximum stays strictly below the 1-ulp ceiling:
    println!("mul F=16: max numerator error = {max_num_err} (< {} ulp fraction)", 1 << 16);
    assert!(max_num_err < 1i128 << 16);
}

#[test]
fn div_floor_bound_f16() {
    for a in EDGES {
        for b in EDGES {
            check_div_f16(a, b);
        }
    }
    let mut rng = Rng(0xB00D_0002);
    for _ in 0..200_000 {
        check_div_f16(rng.i64(), rng.i64());
    }
    // division-specific edges: b = ±1, ±2^16, i64::MIN
    for b in [1i64, -1, 1 << 16, -(1 << 16), i64::MIN] {
        for a in EDGES {
            check_div_f16(a, b);
        }
    }
}

#[test]
fn mul_floor_bound_f32() {
    type F32 = Fixed<32>;
    let mut rng = Rng(0xB00D_0003);
    for _ in 0..50_000 {
        let (a, b) = (rng.i64(), rng.i64());
        let p = (a as i128) * (b as i128);
        match F32::from_raw(a).checked_mul(F32::from_raw(b)) {
            Some(c) => {
                assert_eq!(c.raw() as i128, p >> 32);
                let num_err = p - ((c.raw() as i128) << 32);
                assert!((0..(1i128 << 32)).contains(&num_err));
            }
            None => assert!((p >> 32) < i64::MIN as i128 || (p >> 32) > i64::MAX as i128),
        }
    }
}

#[test]
fn sqrt_bound_f16() {
    // small values exhaustively...
    for a in 0u64..4096 {
        check_sqrt(a as i64);
    }
    let mut rng = Rng(0xB00D_0004);
    // ...then large random samples (non-negative raws)
    for _ in 0..100_000 {
        check_sqrt((rng.next() >> 1) as i64);
    }
    // boundaries around perfect squares
    for k in 0u32..31 {
        let s = 1u64 << ((k + 1) / 2 + 3);
        let sq = s * s;
        if sq < (1 << 62) {
            check_sqrt(sq as i64);
            check_sqrt((sq - 1) as i64);
            check_sqrt((sq + 1) as i64);
        }
    }
}

fn check_sqrt(a: i64) {
    assert!(a >= 0);
    let s_ref = isqrt_ref(a as u64);
    let c = Q16::from_raw(a).sqrt().expect("non-negative, FRAC even");
    assert_eq!(c.raw() as u64, s_ref << 8, "sqrt({a})");
    // Floor property: s_ref^2 <= a < (s_ref+1)^2, hence the worst-case
    // bound of docs/math-spec.md: -2^-8 < err <= 0 for Q16 (<= 256 ulp).
    assert!(
        s_ref * s_ref <= (a as u64) && (a as u64) < (s_ref + 1) * (s_ref + 1),
        "isqrt floor property for {a}"
    );
}

#[test]
fn sqrt_bound_f32() {
    type F32 = Fixed<32>;
    let mut rng = Rng(0xB00D_0005);
    for _ in 0..20_000 {
        let a = (rng.next() >> 1) as i64;
        let s_ref = isqrt_ref(a as u64);
        let c = F32::from_raw(a).sqrt().unwrap();
        assert_eq!(c.raw() as u64, s_ref << 16);
        assert!(s_ref * s_ref <= (a as u64) && (a as u64) < (s_ref + 1) * (s_ref + 1));
    }
}

#[test]
fn pow2_exactness_and_log2_roundtrip_f16() {
    // Exact band: -FRAC <= e <= 62-FRAC => error 0, and log2 is the exact inverse.
    for e in -16i32..=46 {
        let p = Q16::pow2(e).expect("in domain");
        assert_eq!(p.raw(), 1i64 << (16 + e), "e = {e}");
        assert_eq!(p.log2(), Some(e), "log2(pow2({e}))");
    }
    // Below resolution: result 0, error 2^e < 2^-16 (one-sided, documented).
    for e in -40i32..=-17 {
        assert_eq!(Q16::pow2(e), Some(Q16::ZERO), "e = {e}");
    }
    assert_eq!(Q16::pow2(-17).unwrap().log2(), None); // inverse precondition broken
    assert_eq!(Q16::pow2(47), None); // 2^63 does not fit
}

#[test]
fn pow2_log2_f32() {
    type F32 = Fixed<32>;
    for e in -32i32..=30 {
        let p = F32::pow2(e).unwrap();
        assert_eq!(p.raw(), 1i64 << (32 + e), "e = {e}");
        assert_eq!(p.log2(), Some(e));
    }
    assert_eq!(F32::pow2(31), None);
    assert_eq!(F32::pow2(-33), Some(F32::ZERO));
}

#[test]
fn log2_floor_property() {
    let mut rng = Rng(0xB00D_0006);
    for _ in 0..100_000 {
        let raw = rng.i64();
        if raw <= 0 {
            continue;
        }
        let l = Q16::from_raw(raw).log2().unwrap();
        // Exact floor: raw in [2^(l+16), 2^(l+17)), hence -1 < err_log2 <= 0 bits.
        let base = 1u64 << (l + 16);
        let next = 1u64 << (l + 17);
        assert!(
            (base..next).contains(&(raw as u64)),
            "log2({raw}) = {l}: raw outside [2^(l+16), 2^(l+17))"
        );
    }
    // legacy q_log2_i: same contract
    assert_eq!(phase_math::q_log2_i(3), Some(-15));
    assert_eq!(phase_math::q_log2_i(0), None);
    assert_eq!(phase_math::q_log2_i(-5), None);
}

#[test]
fn add_sub_exact_and_invertible() {
    let mut rng = Rng(0xB00D_0007);
    for a in EDGES {
        for b in EDGES {
            check_add_sub(a, b);
        }
    }
    for _ in 0..100_000 {
        check_add_sub(rng.i64(), rng.i64());
    }
}

fn check_add_sub(a: i64, b: i64) {
    // addition: exact (error 0) where defined; sub is the exact inverse
    match Q16::from_raw(a).checked_add(Q16::from_raw(b)) {
        Some(c) => {
            assert_eq!(c.raw(), a.wrapping_add(b), "add {a} + {b}");
            assert_eq!(c.checked_sub(Q16::from_raw(b)), Some(Q16::from_raw(a)));
        }
        None => assert!(a.checked_add(b).is_none(), "add {a} + {b}"),
    }
    // subtraction: exact (error 0) where defined; add is the exact inverse
    match Q16::from_raw(a).checked_sub(Q16::from_raw(b)) {
        Some(c) => {
            assert_eq!(c.raw(), a.wrapping_sub(b), "sub {a} - {b}");
            assert_eq!(c.checked_add(Q16::from_raw(b)), Some(Q16::from_raw(a)));
        }
        None => assert!(a.checked_sub(b).is_none(), "sub {a} - {b}"),
    }
    // negation: MIN is the only hole
    assert_eq!(Q16::from_raw(a).checked_neg().is_none(), a == i64::MIN);
}

#[test]
fn legacy_q_mul_div_bounds_match_docs() {
    let mut rng = Rng(0xB00D_0008);
    // q_mul: conservative domain (None iff raw product does not fit i64),
    // floor rounding: 0 <= p - r*2^16 < 2^16.
    for _ in 0..100_000 {
        let (a, b) = (rng.i64(), rng.i64());
        let p = (a as i128) * (b as i128);
        match phase_math::q_mul(a, b) {
            Some(r) => {
                let num_err = p - ((r as i128) << 16);
                assert!((0..(1i128 << 16)).contains(&num_err), "q_mul {a} {b}");
            }
            None => assert!(p < i64::MIN as i128 || p > i64::MAX as i128, "q_mul {a} {b}"),
        }
        // q_div: truncation toward zero: |num - r*b| < |b|
        if b == 0 {
            assert_eq!(phase_math::q_div(a, b), None);
            continue;
        }
        let num = (a as i128) << 16;
        if let Some(r) = phase_math::q_div(a, b) {
            let rem = num - (r as i128) * (b as i128);
            // truncation: |rem| < |b| and rem has the sign of num (or 0)
            assert!(rem.abs() < (b as i128).abs(), "q_div {a} {b}");
            assert!(rem == 0 || (rem > 0) == (num > 0), "q_div sign {a} {b}");
        }
    }
    // conservative domain demo: product overflows i64 although the shifted
    // result would fit — documented behavior of the legacy check
    assert_eq!(phase_math::q_mul(1 << 40, 1 << 24), None);
    // ((1<<40)*(1<<24)) >> 16 = 2^48 would fit; the raw product 2^64 does not
}
