//! Bijections on Z/2^64 with analytically computed inverses.
//!
//! Every operation in this module is an exact bijection (error 0): applying
//! the analytic inverse after the forward op returns the original word for
//! every input. Composition of bijections is a bijection, so a `Compose`
//! pipeline inverts by reversing and negating its ops.

/// Modular inverse of an odd `a` modulo 2^64 (Newton iteration, 6 steps).
///
/// Domain: `a` must be odd (debug-asserted). Exact: `a * mod_inv_odd(a) ≡ 1
/// (mod 2^64)`. Convergence: `x0 = a` is already correct mod 8 (3 bits, odd
/// `a` squares to 1 mod 8) and Newton doubles the correct bits each step,
/// so 6 iterations give 3·2^6 = 192 ≥ 64 bits. The map `a ↦ a⁻¹` is itself
/// a bijection on the odd residues and an involution.
pub fn mod_inv_odd(a: u64) -> u64 {
    debug_assert!(a & 1 == 1);
    let mut x = a;
    for _ in 0..6 {
        x = x.wrapping_mul(2u64.wrapping_sub(a.wrapping_mul(x)));
    }
    x
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    /// `x = x * m`, `m` odd (bijection; inverse is `MulOdd(mod_inv_odd(m))`)
    MulOdd(u64),
    /// `x = x + c` (wrapping addition; inverse is `Add(-c)`)
    Add(u64),
    /// `x = x ^ k` (inverse is itself: an involution)
    Xor(u64),
    /// `x = rotl(x, r)` (`r` is taken mod 64; inverse is `RotL(64 - r)`)
    RotL(u32),
}

impl Op {
    /// Applies the operation to `x`. Exact bijection on Z/2^64 (error 0).
    pub fn apply(&self, x: u64) -> u64 {
        match *self {
            Op::MulOdd(m) => x.wrapping_mul(m),
            Op::Add(c) => x.wrapping_add(c),
            Op::Xor(k) => x ^ k,
            Op::RotL(r) => x.rotate_left(r & 63),
        }
    }

    /// The analytic inverse operation: `self.inverse().apply(self.apply(x)) == x`
    /// for every `x`. Exact (error 0).
    pub fn inverse(&self) -> Op {
        match *self {
            Op::MulOdd(m) => Op::MulOdd(mod_inv_odd(m)),
            Op::Add(c) => Op::Add(c.wrapping_neg()),
            Op::Xor(k) => Op::Xor(k),
            Op::RotL(r) => Op::RotL((64 - (r & 63)) & 63),
        }
    }
}

/// A sequence of bijections: the composition is itself a bijection and the
/// inverse is analytic (reverse the ops, invert each).
#[derive(Debug, Clone, Default)]
pub struct Compose {
    ops: Vec<Op>,
}

impl Compose {
    pub fn new() -> Self {
        Compose { ops: Vec::new() }
    }

    /// Appends one bijective op to the pipeline.
    pub fn push(&mut self, op: Op) {
        self.ops.push(op);
    }

    /// Applies the pipeline to `x` in push order. Exact bijection (error 0).
    pub fn apply(&self, mut x: u64) -> u64 {
        for o in &self.ops {
            x = o.apply(x);
        }
        x
    }

    /// The exact inverse pipeline: `self.inverse().apply(self.apply(x)) == x`
    /// for every `x` (composition of bijections, reversed and negated).
    pub fn inverse(&self) -> Compose {
        Compose {
            ops: self.ops.iter().rev().map(|o| o.inverse()).collect(),
        }
    }

    pub fn len(&self) -> usize {
        self.ops.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ops.is_empty()
    }
}

/// Reversible ARX round on a word pair (bijection of (u64, u64)).
///
/// Exact: `arx_round_inv(arx_round(a, b, k), k) == (a, b)` for all inputs;
/// `k` is taken mod 64 on both sides.
pub fn arx_round(a: u64, b: u64, k: u32) -> (u64, u64) {
    let s = a.wrapping_add(b);
    let t = b.rotate_left(k & 63) ^ s;
    (s, t)
}

/// Exact inverse of [`arx_round`].
pub fn arx_round_inv(s: u64, t: u64, k: u32) -> (u64, u64) {
    let b = (t ^ s).rotate_right(k & 63);
    let a = s.wrapping_sub(b);
    (a, b)
}
