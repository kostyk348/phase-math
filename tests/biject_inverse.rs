//! Inverse-property tests for the exact bijection family
//! (docs/math-spec.md, Family A): `inverse(forward(x)) == x` must hold for
//! every input — verified here on edge values and large random samples.

use phase_math::{arx_round, arx_round_inv, mod_inv_odd, Compose, Op};

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
}

const EDGES: [u64; 8] = [
    0,
    1,
    2,
    u64::MAX,
    u64::MAX - 1,
    1 << 63,
    (1 << 63) + 1,
    0xDEAD_BEEF_CAFE_F00D,
];

#[test]
fn mul_odd_inverse_and_multiplier_identity() {
    let mut rng = Rng(0xB10C_0301);
    for i in 0..50_000u64 {
        let m = (if i == 0 { 3 } else { rng.next() }) | 1; // odd
        let op = Op::MulOdd(m);
        // multiplier identity: m * m^-1 ≡ 1 (mod 2^64)
        assert_eq!(m.wrapping_mul(mod_inv_odd(m)), 1, "m = {m:#x}");
        for x in EDGES {
            assert_eq!(op.inverse().apply(op.apply(x)), x, "m={m:#x} x={x:#x}");
        }
        for _ in 0..20 {
            let x = rng.next();
            assert_eq!(op.inverse().apply(op.apply(x)), x);
        }
    }
}

#[test]
fn add_xor_rotl_inverse() {
    let mut rng = Rng(0xB10C_0302);
    for i in 0..50_000u64 {
        let c = if i == 0 { u64::MAX } else { rng.next() };
        let k = if i == 1 { u64::MAX } else { rng.next() };
        let r = (rng.next() % 64) as u32;
        for x in EDGES {
            assert_eq!(Op::Add(c).inverse().apply(Op::Add(c).apply(x)), x);
            assert_eq!(Op::Xor(k).inverse().apply(Op::Xor(k).apply(x)), x);
            assert_eq!(Op::RotL(r).inverse().apply(Op::RotL(r).apply(x)), x);
        }
    }
    // every rotation amount, exhaustively on edges
    for r in 0u32..64 {
        for x in EDGES {
            assert_eq!(Op::RotL(r).inverse().apply(Op::RotL(r).apply(x)), x);
        }
    }
}

#[test]
fn mod_inv_odd_is_involution() {
    let mut rng = Rng(0xB10C_0303);
    // the exact inverse of mod_inv_odd is itself: (a^-1)^-1 == a
    assert_eq!(mod_inv_odd(mod_inv_odd(1)), 1);
    assert_eq!(mod_inv_odd(mod_inv_odd(u64::MAX)), u64::MAX);
    for _ in 0..50_000 {
        let a = rng.next() | 1;
        assert_eq!(mod_inv_odd(mod_inv_odd(a)), a);
    }
}

#[test]
fn compose_inverse_random_pipelines() {
    let mut rng = Rng(0xB10C_0304);
    for _ in 0..5_000 {
        let mut c = Compose::new();
        let n = 1 + rng.next() % 8;
        for _ in 0..n {
            let op = match rng.next() % 4 {
                0 => Op::MulOdd(rng.next() | 1),
                1 => Op::Add(rng.next()),
                2 => Op::Xor(rng.next()),
                _ => Op::RotL((rng.next() % 64) as u32),
            };
            c.push(op);
        }
        let x = rng.next();
        assert_eq!(c.inverse().apply(c.apply(x)), x, "pipeline len {}", c.len());
        assert_eq!(c.is_empty(), false);
    }
    // empty pipeline is the identity bijection
    let e = Compose::new();
    assert_eq!(e.inverse().apply(e.apply(7)), 7);
}

#[test]
fn arx_round_inverse() {
    let mut rng = Rng(0xB10C_0305);
    // edges: (0,0), (MAX, MAX), k = 0 and 63
    for k in [0u32, 1, 31, 63, 64, 65, 127] {
        for (a, b) in [(0u64, 0u64), (u64::MAX, u64::MAX), (0, u64::MAX), (1, 2)] {
            assert_eq!(
                {
                    let (s, t) = arx_round(a, b, k);
                    arx_round_inv(s, t, k)
                },
                (a, b),
                "k={k} a={a} b={b}"
            );
        }
    }
    for _ in 0..100_000 {
        let a = rng.next();
        let b = rng.next();
        let k = (rng.next() % 128) as u32; // includes k >= 64: masked mod 64
        let (s, t) = arx_round(a, b, k);
        assert_eq!(arx_round_inv(s, t, k), (a, b));
    }
}
