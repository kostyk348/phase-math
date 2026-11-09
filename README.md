# phase-math — reversible math for the PHASE paradigm

English | [Русский](README.ru.md)

Math **inside a phase** = bijections on Z/2^64 with **analytically computed
inverses**: `MulOdd(m)` (multiplying by an odd constant is a bijection; the
inverse is `MulOdd(m⁻¹)` via Newton mod-inverse), `Add`, `Xor`, `RotL`,
reversible ARX rounds. Any composition is a bijection, so any formula is
exactly undoable — `F⁻¹(F(x)) == x` is a property, not a promise.

Float math is not invertible and stays **at the phase boundary**. This crate
is the numerical layer for that boundary: fixed-point and interval arithmetic
with documented, test-verified error budgets, so a VM can decide which
computations can be reversed exactly and what a rollback must undo.

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![no deps](https://img.shields.io/badge/dependencies-0-green.svg)]()
[![Release](https://img.shields.io/github/v/release/kostyk348/phase-math?sort=semver&color=blue)](https://github.com/kostyk348/phase-math/releases)
[![Last commit](https://img.shields.io/github/last-commit/kostyk348/phase-math?branch=main)](https://github.com/kostyk348/phase-math/commits/main)

## Quickstart

phase-math is not published on crates.io yet — add it from git:

```bash
cargo add --git https://github.com/kostyk348/phase-math phase-math
```

or in `Cargo.toml`:

```toml
[dependencies]
phase-math = { git = "https://github.com/kostyk348/phase-math" }
```

Minimal example:

```rust
use phase_math::{Compose, Interval, Op, Q16};

fn main() {
    // 1. Bijections on Z/2^64: any composition inverts analytically.
    let mut enc = Compose::new();
    enc.push(Op::MulOdd(0x9E37_79B9_7F4A_7C15)); // odd multiplier
    enc.push(Op::Add(0xD1B5_4A32_D192_ED03));
    let dec = enc.inverse();
    assert_eq!(dec.apply(enc.apply(42)), 42);

    // 2. Fixed point (Q16.16 here): exact when representable,
    //    None (or a documented error bound) otherwise.
    let x = Q16::from_i32(3).unwrap(); // 3.0
    let y = Q16::from_i32(7).unwrap();
    assert_eq!(x.checked_mul(y), Q16::from_i32(21)); // exact: 21 * 2^16 divides
    assert_eq!(x.checked_mul(Q16::MAX), None); // out of domain -> roll back

    // 3. Saturating variants report clamping with a flag.
    let s = Q16::MAX.saturating_add(Q16::ONE);
    assert!(s.saturated && s.value == Q16::MAX);

    // 4. Intervals with outward rounding: the exact result stays inside.
    let around_one = Interval::new(Q16::from_raw(65535), Q16::from_raw(65537)).unwrap();
    let p = around_one.mul(Interval::point(x)).unwrap(); // 3.0 * (1 ± 1 ulp)
    assert!(p.contains(x));
    println!("3.0 * [1 ∓ 1 ulp] = [{}, {}]", p.lo, p.hi);
    // -> [2.9999542236328125, 3.0000457763671875]
}
```

## What is inside

- `Op::{MulOdd, Add, Xor, RotL}` — bijective primitives with
  `.apply()` / `.inverse()`; `mod_inv_odd(a)` — modular inverse mod 2^64;
  `Compose` — builder with exact analytic inversion;
  `arx_round` / `arx_round_inv` — invertible ARX round on a word pair;
- `Fixed<FRAC>` — generic fixed point (value = `raw / 2^FRAC`, `i64` raw);
  `Q16` is the classic Q16.16. Checked ops return `Option` (an out-of-domain
  result is `None` and the caller rolls the phase back — fail = un-run);
  saturating ops return the value plus a `saturated` flag;
  `sqrt`, `pow2`, `log2`, exact decimal `Display`;
- `Interval<FRAC>` — closed intervals with outward rounding: the exact real
  result is guaranteed to stay inside, or the operation returns `None`;
- legacy `q_*` free functions over `Q = i64` (Q16.16), kept for compatibility.

## Error budgets and classification

Every operation is classified and budgeted in
[docs/math-spec.md](docs/math-spec.md): operation → domain → reverse
operation → class (bijective / partially invertible / lossy) → error bound
or containment guarantee. Examples:

| Operation | Class | Worst-case error |
|---|---|---|
| `MulOdd`, `Add`, `Xor`, `RotL`, `Compose`, `arx_round` | bijection | exact (0) |
| `Fixed::checked_mul` | lossy (floor) | `−2^−FRAC < err ≤ 0` (≤ 1 ulp) |
| `Fixed::checked_div` | lossy (floor) | `−2^−FRAC < err ≤ 0` (≤ 1 ulp) |
| `Fixed::sqrt` | lossy (floor) | `−2^−(FRAC/2) < err ≤ 0` (≤ 256 ulp for Q16) |
| `Interval::mul` / `div` | outward per endpoint | endpoints within 1 ulp of the exact extremes; width ≤ exact + 2 ulp |
| saturating ops | lossy when clamped | flag `true` = no numeric bound, value is the clamp |

## Honesty and scope

- Reversible-friendly arithmetic for VM use: the bijective family is exactly
  invertible; the fixed-point and interval families are domain-checked and
  carry error budgets so a VM knows what a rollback must undo.
- The error bounds are **worst-case upper bounds**, derived analytically and
  enforced by deterministic randomized tests against exact i128 references
  (no floats anywhere, including the tests). They are test-verified
  guarantees, not machine-checked formal proofs.
- No promises beyond implemented features: no floats, no transcendentals
  beyond `2^x` and `floor(log2)`, no NaN or infinities; intervals return
  `None` instead of inflating; saturating results carry a flag instead of a
  bound; the legacy `q_*` API is documented exactly as it behaves
  (including its conservative overflow check) and is frozen.
- The API may still evolve before 1.0. The specification in
  [docs/math-spec.md](docs/math-spec.md) is versioned with the crate.

## Tests

56 tests total (`cargo test`): bijective inverses on edge values and large
random samples; error-bound tests against exact i128 references (200k random
full-range samples for mul/div bounds); saturating budget checks; interval
containment and exact-outward-endpoint tightness; partial-inverse
preconditions and concrete lossy examples; doc tests for the documented
examples. The full claim-to-test map is in
[docs/math-spec.md](docs/math-spec.md).

Part of the PHASE architecture (github.com/kostyk348/phase-arch) — the
reversible-math layer; reference implementation lives with phase-vm.
