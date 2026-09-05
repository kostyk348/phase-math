//! Биекции на Z/2^64 с аналитическими инверсиями.

/// Обратное к нечётному a по модулю 2^64 (Newton: 6 итераций).
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
    /// x = x * m, m нечётное (биекция)
    MulOdd(u64),
    /// x = x + c (wrap)
    Add(u64),
    /// x = x ^ k
    Xor(u64),
    /// x = rotl(x, r)
    RotL(u32),
}

impl Op {
    pub fn apply(&self, x: u64) -> u64 {
        match *self {
            Op::MulOdd(m) => x.wrapping_mul(m),
            Op::Add(c) => x.wrapping_add(c),
            Op::Xor(k) => x ^ k,
            Op::RotL(r) => x.rotate_left(r & 63),
        }
    }
    pub fn inverse(&self) -> Op {
        match *self {
            Op::MulOdd(m) => Op::MulOdd(mod_inv_odd(m)),
            Op::Add(c) => Op::Add(c.wrapping_neg()),
            Op::Xor(k) => Op::Xor(k),
            Op::RotL(r) => Op::RotL((64 - (r & 63)) & 63),
        }
    }
}

/// Последовательность биекций: композиция = биекция, инверсия аналитическая.
#[derive(Debug, Clone, Default)]
pub struct Compose {
    ops: Vec<Op>,
}

impl Compose {
    pub fn new() -> Self {
        Compose { ops: Vec::new() }
    }
    pub fn push(&mut self, op: Op) {
        self.ops.push(op);
    }
    pub fn apply(&self, mut x: u64) -> u64 {
        for o in &self.ops {
            x = o.apply(x);
        }
        x
    }
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

/// Обратимый ARX-раунд (биекция пары слов).
pub fn arx_round(a: u64, b: u64, k: u32) -> (u64, u64) {
    let s = a.wrapping_add(b);
    let t = b.rotate_left(k & 63) ^ s;
    (s, t)
}
pub fn arx_round_inv(s: u64, t: u64, k: u32) -> (u64, u64) {
    let b = (t ^ s).rotate_right(k & 63);
    let a = s.wrapping_sub(b);
    (a, b)
}
