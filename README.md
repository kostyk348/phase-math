# phase-math — reversible math for the PHASE paradigm

Math **inside a phase** = bijections on Z/2^64 with **analytically computed
inverses**: `MulOdd(m)` (multiplying by an odd constant is a bijection; the
inverse is `MulOdd(m⁻¹)` via Newton mod-inverse), `Add`, `Xor`, `RotL`,
reversible ARX rounds. Any composition is a bijection, so any formula is
exactly undoable — `F⁻¹(F(x)) == x` is a property, not a promise.

Float math is not invertible and stays **at the phase boundary**.

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![no deps](https://img.shields.io/badge/dependencies-0-green.svg)]()
[![Release](https://img.shields.io/github/v/release/kostyk348/phase-math?sort=semver&color=blue)](https://github.com/kostyk348/phase-math/releases)
[![Last commit](https://img.shields.io/github/last-commit/kostyk348/phase-math?branch=main)](https://github.com/kostyk348/phase-math/commits/main)

## API

- `Op::{MulOdd,Add,Xor,RotL}` — bijective primitives with `.apply()`/`.inverse()`;
- `mod_inv_odd(a)` — modular inverse mod 2^64 (Newton, 6 iterations);
- `Compose` — builder; `.apply(x)` / `.inverse()` (exact);
- `arx_round` / `arx_round_inv` — invertible ARX round on a word pair;
- Q16.16 (`Q`): `q_from_i32`, `q_mul`, `q_div`, `q_sqrt` — fixed-point with
  **boundary domain validation**: they return `Option`, so an out-of-domain
  result is `None` and the caller rolls the phase back (fail = un-run).

```rust
let mut enc = Compose::new();
enc.push(Op::MulOdd(0x9E37_79B9_7F4A_7C15)); // odd multiplier
enc.push(Op::Add(0xD1B5_4A32_D192_ED03));
let dec = enc.inverse();
assert_eq!(dec.apply(enc.apply(x)), x);
```

## Tests

mod-inverse, compose round-trip (3000 random), ARX invertibility,
Q16.16 mul/div/sqrt precision + domain checks.

Part of the PHASE architecture (github.com/kostyk348/phase-arch) — the
reversible-math layer; reference implementation lives with phase-vm.
