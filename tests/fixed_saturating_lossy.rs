//! Saturating budget checks, lossy-operation demonstrations and
//! partial-inverse preconditions (docs/math-spec.md, Families B, C, D).

use phase_math::{q_div, q_log2_i, q_mul, q_pow2, q_sqrt, Q16};

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

const EDGES: [i64; 8] = [0, 1, -1, i64::MAX, i64::MIN, i64::MAX - 1, i64::MIN + 1, 1 << 16];

#[test]
fn saturating_add_sub_budget() {
    let mut rng = Rng(0x5A7_0001);
    // Contract: flag == false => exact result (error 0); flag == true =>
    // clamped to MAX/MIN (no numeric bound; the flag is the signal).
    for _ in 0..100_000 {
        let (a, b) = (rng.i64(), rng.i64());
        let s = Q16::from_raw(a).saturating_add(Q16::from_raw(b));
        match a.checked_add(b) {
            Some(exact) => {
                assert!(!s.saturated, "add {a} {b}");
                assert_eq!(s.value.raw(), exact);
            }
            None => {
                assert!(s.saturated, "add {a} {b}");
                assert!(s.value.raw() == i64::MAX || s.value.raw() == i64::MIN);
                assert_eq!(s.value.raw(), a.saturating_add(b));
            }
        }
        let d = Q16::from_raw(a).saturating_sub(Q16::from_raw(b));
        match a.checked_sub(b) {
            Some(exact) => {
                assert!(!d.saturated, "sub {a} {b}");
                assert_eq!(d.value.raw(), exact);
            }
            None => {
                assert!(d.saturated, "sub {a} {b}");
                assert!(d.value.raw() == i64::MAX || d.value.raw() == i64::MIN);
                assert_eq!(d.value.raw(), a.saturating_sub(b));
            }
        }
    }
}

#[test]
fn saturating_mul_budget_vs_i128_reference() {
    let mut rng = Rng(0x5A7_0002);
    for a in EDGES {
        for b in EDGES {
            check_sat_mul(a, b);
        }
    }
    for _ in 0..100_000 {
        check_sat_mul(rng.i64(), rng.i64());
    }
}

fn check_sat_mul(a: i64, b: i64) {
    let s = Q16::from_raw(a).saturating_mul(Q16::from_raw(b));
    let r = ((a as i128) * (b as i128)) >> 16; // exact floor reference
    let clamped = r.clamp(i64::MIN as i128, i64::MAX as i128);
    assert_eq!(s.value.raw() as i128, clamped, "sat_mul {a} {b}");
    assert_eq!(s.saturated, clamped != r, "sat_mul flag {a} {b}");
    if !s.saturated {
        // unsaturated budget: exactly the checked-mul result
        assert_eq!(s.value.raw() as i128, r);
        assert_eq!(
            s.value,
            Q16::from_raw(a).checked_mul(Q16::from_raw(b)).unwrap()
        );
    }
}

#[test]
fn saturating_lossy_concrete_pairs() {
    // Two different addends clamp to the same value: the addend is
    // unrecoverable. f(x) = MAX.saturating_add(x):
    assert_eq!(Q16::MAX.saturating_add(Q16::ONE).value, Q16::MAX);
    assert_eq!(Q16::MAX.saturating_add(Q16::from_raw(2)).value, Q16::MAX);
    // Symmetric at the floor: g(x) = MIN.saturating_sub(x):
    assert_eq!(Q16::MIN.saturating_sub(Q16::ONE).value, Q16::MIN);
    assert_eq!(Q16::MIN.saturating_sub(Q16::from_raw(2)).value, Q16::MIN);
    // Saturating mul collapses two multipliers at the ceiling:
    assert_eq!(
        Q16::MAX.saturating_mul(Q16::from_raw(1 << 18)).value,
        Q16::MAX
    );
    assert_eq!(
        Q16::MAX.saturating_mul(Q16::from_raw(1 << 19)).value,
        Q16::MAX
    );
}

#[test]
fn mul_floor_lossy_concrete_pairs() {
    // Fixed mul (floor): (3 * 2^-16) * (1 * 2^-16) and (3 * 2^-16) * (2 * 2^-16)
    // both floor to 0 — the second factor is destroyed:
    assert_eq!(
        Q16::from_raw(3).checked_mul(Q16::from_raw(1)),
        Some(Q16::ZERO)
    );
    assert_eq!(
        Q16::from_raw(3).checked_mul(Q16::from_raw(2)),
        Some(Q16::ZERO)
    );
    // legacy q_mul, same pair:
    assert_eq!(q_mul(3, 1), Some(0));
    assert_eq!(q_mul(3, 2), Some(0));
    // division: two numerators collapse (legacy truncs toward 0, Fixed floors):
    assert_eq!(q_div(1, 1 << 20), Some(0));
    assert_eq!(q_div(2, 1 << 20), Some(0));
    assert_eq!(
        Q16::from_raw(1).checked_div(Q16::from_raw(1 << 20)),
        Some(Q16::ZERO)
    );
    assert_eq!(
        Q16::from_raw(2).checked_div(Q16::from_raw(1 << 20)),
        Some(Q16::ZERO)
    );
}

#[test]
fn sqrt_log2_lossy_concrete_pairs() {
    // sqrt: 16.0 and the next representable value above it both map to 4.0
    assert_eq!(Q16::from_raw(1 << 20).sqrt(), Some(Q16::from_raw(1 << 18)));
    assert_eq!(
        Q16::from_raw((1 << 20) + 1).sqrt(),
        Some(Q16::from_raw(1 << 18))
    );
    assert_eq!(q_sqrt(1 << 20), Some(1 << 18));
    assert_eq!(q_sqrt((1 << 20) + 1), Some(1 << 18));
    // log2: 1.0 and 1.0 + 1 ulp both map to 0
    assert_eq!(Q16::from_raw(1 << 16).log2(), Some(0));
    assert_eq!(Q16::from_raw((1 << 16) + 1).log2(), Some(0));
    assert_eq!(q_log2_i((1 << 16) + 1), Some(0));
}

#[test]
fn partial_inverse_mul_div_precondition() {
    // Documented precondition (Family B): q_div(q_mul(a, b), b) == a exactly
    // IFF the raw product a*b is divisible by 2^16 (and fits in i64, b != 0).
    let mut rng = Rng(0x5A7_0003);
    for _ in 0..50_000 {
        let a = (rng.next() % (1 << 30)) as i64 - (1 << 29);
        let k = (rng.next() % (1 << 10)) as i64 - (1 << 9);
        if k == 0 {
            continue;
        }
        let b = k << 16; // multiple of 2^16 => product divisible by 2^16
        let p = q_mul(a, b).expect("in domain");
        assert_eq!(p, a * k); // mul exact under the precondition (Q units)
        assert_eq!(q_div(p, b), Some(a), "exact inverse under precondition");
    }
    // Counterexample: without divisibility the roundtrip is lossy.
    let (a, b) = (3i64, 1i64);
    let p = q_mul(a, b).unwrap();
    assert_eq!(p, 0);
    assert_ne!(q_div(p, b), Some(a)); // 0/1 = 0 != 3: information lost
}

#[test]
fn partial_inverse_pow2_log2_precondition() {
    // Legacy pair: exact inverse exactly on e in [-16, 31].
    for e in -16i32..=31 {
        assert_eq!(q_log2_i(q_pow2(e).unwrap()), Some(e), "e = {e}");
    }
    // Outside: pow2 collapses (2^-17 < 1 ulp) and the inverse is undefined.
    assert_eq!(q_pow2(-17), Some(0));
    assert_eq!(q_log2_i(0), None);
    assert!(q_pow2(-32).is_none());
    assert!(q_pow2(32).is_none());
    // Fixed pair: exact inverse on e in [-FRAC, 62-FRAC].
    for e in -16i32..=46 {
        assert_eq!(Q16::pow2(e).unwrap().log2(), Some(e), "e = {e}");
    }
    assert_eq!(Q16::pow2(-17).unwrap().log2(), None);
}

#[test]
fn fixed_add_sub_partial_inverse_documented() {
    // Precondition: every intermediate stays in i64 (checked ops return None
    // otherwise, and the caller rolls the phase back).
    let mut rng = Rng(0x5A7_0004);
    for _ in 0..50_000 {
        let (a, b) = (rng.i64(), rng.i64());
        if let Some(c) = Q16::from_raw(a).checked_add(Q16::from_raw(b)) {
            assert_eq!(c.checked_sub(Q16::from_raw(b)), Some(Q16::from_raw(a)));
        }
        if let Some(c) = Q16::from_raw(a).checked_sub(Q16::from_raw(b)) {
            assert_eq!(c.checked_add(Q16::from_raw(b)), Some(Q16::from_raw(a)));
        }
    }
}

#[test]
fn from_i32_inverse_precondition() {
    // from_i32 is exact and injective; the inverse (raw >> 16) is exact iff
    // the low FRAC raw bits are zero.
    for v in [0i32, 1, -1, i32::MAX, i32::MIN, 12345] {
        let x = Q16::from_i32(v).unwrap();
        assert_eq!(x.raw() >> 16, v as i64); // exact inverse on the image
        assert_eq!(x.raw() & 0xFFFF, 0); // low bits zero by construction
    }
    // a raw with nonzero low bits does not come back through >> 16:
    assert_ne!(Q16::from_raw(1).raw() >> 16, 1);
}
