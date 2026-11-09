# phase-math — mathematical specification

This document specifies every arithmetic operation of phase-math (matching
crate version 0.2.0): its domain, its reverse operation if one exists, its
invertibility class, and its worst-case error bound or containment
guarantee. Every row was cross-checked against the source of this commit
and is enforced by the test named in
[How the bounds are tested](#how-the-bounds-are-tested).

The classification vocabulary is the one phase-vm consumes: it decides
which computations can be reversed exactly (bijective), under which
precondition (partially invertible), and where information is destroyed
(lossy) so a phase boundary is required.

## Notation and conventions

- A fixed-point value with `FRAC` fractional bits is `x = raw · 2^-FRAC`
  with `raw ∈ i64`. The legacy family fixes `FRAC = 16` (`type Q = i64`,
  the classic Q16.16 layout: 48 integer bits including sign, 16 fraction
  bits). One **ulp** (unit in the last place) is `2^-FRAC`.
- `err` denotes `computed − exact real result`.
- A bound `err ∈ (−1 ulp, 0]` means the operation never overestimates and
  misses the exact value by strictly less than one ulp (one-sided,
  rounding toward −infinity). Two-sided bounds are written `|err| < 1 ulp`.
- All bounds below are **derived** (see derivations inline) and **verified
  by tests** against exact i128 rational references. There are no floats
  anywhere in the crate or its tests.

## Invertibility classes

- **bijection** — the analytic inverse is provided and is exact over the
  whole stated domain; error is 0.
- **partially invertible** — an exact inverse exists on a documented
  subdomain (the precondition); outside it the operation either refuses
  (`None` / flag) or collapses inputs. The precondition is part of the
  contract.
- **lossy** — no inverse exists (many-to-one). A worst-case error bound
  and/or a concrete example pair demonstrating the information loss is
  given instead.

## Family A — exact bijections on Z/2^64 (`src/biject.rs`)

All operations are total, exact (`err = 0`) and carry analytic inverses;
composition of any of them is again a bijection.

| Operation | Domain | Reverse operation | Class | Error |
|---|---|---|---|---|
| `Op::MulOdd(m).apply(x)` | all `x ∈ u64`; `m` odd | `Op::MulOdd(mod_inv_odd(m))` | bijection | 0 |
| `Op::Add(c).apply(x)` | all `x ∈ u64` | `Op::Add(c.wrapping_neg())` | bijection | 0 |
| `Op::Xor(k).apply(x)` | all `x ∈ u64` | `Op::Xor(k)` (involution) | bijection | 0 |
| `Op::RotL(r).apply(x)` | all `x ∈ u64`; `r` used mod 64 | `Op::RotL((64 − (r mod 64)) mod 64)` | bijection | 0 |
| `mod_inv_odd(a)` | `a` odd | itself (`(a⁻¹)⁻¹ = a`) | bijection on odd residues | 0 |
| `Compose::apply` (any pipeline) | all `x ∈ u64` | `Compose::inverse()` — ops reversed and inverted | bijection | 0 |
| `arx_round(a, b, k)` | all `(u64, u64)`; `k` used mod 64 | `arx_round_inv(s, t, k)` | bijection on word pairs | 0 |

Why these are bijections:

- `x ↦ x·m mod 2^64` is injective iff `gcd(m, 2^64) = 1` iff `m` is odd.
  `mod_inv_odd` computes the inverse with 6 Newton steps: `x0 = a` is
  already correct mod 8 (any odd `a` satisfies `a·a ≡ 1 mod 8`), Newton
  doubles the correct bits per step, so 6 steps give `3·2^6 = 192 ≥ 64`
  bits; the identity `a·a⁻¹ ≡ 1 (mod 2^64)` is test-enforced.
- `arx_round`: `(a, b) ↦ (s, t) = (a + b, rotl(b, k) ^ s)`; the inverse
  recomputes `b = rotr(t ^ s, k)` then `a = s − b` — subtraction undoes
  the addition, rotate-xor undoes itself.

Tests: `tests/biject_inverse.rs` — per-operation
`inverse(forward(x)) == x` on edge words (0, 1, `u64::MAX`, …) and random
samples (50k ops with 20 random inputs each, all 64 rotation amounts,
100k ARX pairs including `k ≥ 64`, 5k random pipelines up to 8 ops), plus
the multiplier identity and the `mod_inv_odd` involution.

## Family B — legacy Q16.16 free functions (`src/q.rs`; `Q = i64`, 1 ulp = 2^−16)

| Operation | Domain (`None` outside) | Reverse operation | Class | Error bound |
|---|---|---|---|---|
| `q_from_i32(v)` | all `i32` | `raw >> 16`, exact on the image (low 16 raw bits zero) | bijection onto image | 0 |
| `q_mul(a, b)` | raw product `a·b` fits `i64` — conservative, see note | — | lossy (floor) | `−2^−16 < err ≤ 0` |
| `q_div(a, b)` | `b ≠ 0` and `(a·2^16)/b` fits `i64` | — | lossy (truncation toward 0) | `|err| < 2^−16` |
| `q_sqrt(a)` | `a ≥ 0` | — | lossy (floor) | `−2^−8 < err ≤ 0` |
| `q_pow2(e)` | `e ∈ [−31, 31]` | `q_log2_i`, exact on `e ∈ [−16, 31]` | partially invertible | 0 for `e ≥ −16`; `err = 2^e < 2^−16` for `e ≤ −17` |
| `q_log2_i(x)` | `x > 0` | — | lossy (exact floor) | `−1 < err_log2 ≤ 0` (bits) |

Notes:

- **Conservative domain of `q_mul`**: it returns `None` whenever the raw
  product `a·b` itself does not fit in `i64`, even when the shifted result
  would (example kept as a test: `q_mul(2^40, 2^24)` is `None`, although
  `2^64 >> 16 = 2^48` is representable). This is preserved for backward
  compatibility; `Fixed::checked_mul` implements the tight rule instead.
- **`q_div` rounds toward zero** (Rust integer division), so its error is
  two-sided. `Fixed::checked_div` floors; the difference is deliberate and
  both behaviors are documented and tested.
- **`q_sqrt` loses half the raw precision** (squaring halves the
  exponent): the bound is 256 ulp in absolute terms for Q16.16.

## Family C — `Fixed<FRAC>` checked operations (`src/fixed.rs`)

| Operation | Domain (`None` outside) | Reverse operation | Class | Error bound |
|---|---|---|---|---|
| `checked_add(a, b)` | no `i64` raw overflow | `a.checked_add(b)?.checked_sub(b) == a` | partially invertible (precondition: no overflow) | 0 |
| `checked_sub(a, b)` | no `i64` raw overflow | `checked_add` (exact) | partially invertible | 0 |
| `checked_neg` | `raw ≠ i64::MIN` | itself | partially invertible | 0 |
| `from_i32(v)` | `v·2^FRAC` fits `i64` | `raw >> FRAC`, exact on the image | bijection onto image | 0 |
| `checked_mul` | `floor(a·b / 2^FRAC)` fits `i64` (tight rule) | — | lossy (floor) | `−2^−FRAC < err ≤ 0` |
| `checked_div` | `b ≠ 0`, `floor(a·2^FRAC / b)` fits `i64` | — | lossy (floor, any sign of `b`) | `−2^−FRAC < err ≤ 0` |
| `sqrt` | `raw ≥ 0`, `FRAC` even | — | lossy (floor) | `−2^−(FRAC/2) < err ≤ 0` |
| `pow2(e)` | `e ≤ 62 − FRAC` (returns 0 for `e < −FRAC`) | `log2`, exact on `e ∈ [−FRAC, 62 − FRAC]` | partially invertible | 0 for `e ≥ −FRAC`; `err = 2^e < 2^−FRAC` for `e < −FRAC` |
| `log2` | `raw > 0` | — | lossy (exact floor) | `−1 < err_log2 ≤ 0` (bits) |

Derivations:

- **mul**: exact = `a·b·2^−2F`, computed = `floor(a·b/2^F)·2^−F`, so
  `err = (floor(a·b/2^F) − a·b/2^F)·2^−F ∈ (−2^−F, 0]` (arithmetic shift
  on i128 is floor division). Error is 0 exactly when `2^F` divides
  `a·b`.
- **div**: `err = (floor(a·2^F/b) − a·2^F/b)·2^−F`. The floor-division
  remainder has the sign of the divisor; dividing by a negative `b`
  flips it back, so the bound is one-sided for **every** `b ≠ 0`
  (verified with negative divisors in the tests).
- **sqrt**: `err = (isqrt(a) − sqrt(a))·2^−(F/2)`; `isqrt` is the true
  floor (`isqrt(a)² ≤ a < (isqrt(a)+1)²`, test-enforced), hence
  `sqrt(a) − isqrt(a) ∈ [0, 1)`.
- **pow2**: exact while the raw `2^(F+e) ≥ 1`; below resolution the exact
  value is under 1 ulp and the floor is 0 with `err = 2^e < 2^−F`.

Tests: `tests/fixed_bounds.rs` — `mul_floor_bound_f16` (200k random
full-range pairs + edge matrix, plus a 50k-sample `FRAC = 32` variant),
`div_floor_bound_f16` (200k + edge divisors ±1, ±2^16, `i64::MIN`),
`sqrt_bound_f16` (100k random non-negative raws + exhaustive small values
+ perfect-square boundaries, against an independent binary-search
reference), `pow2_exactness_and_log2_roundtrip_f16`, `log2_floor_property`
(100k samples), `add_sub_exact_and_invertible` (100k + edges, both
directions of the inverse), and `legacy_q_mul_div_bounds_match_docs` for
Family B.

## Family D — `Fixed<FRAC>` saturating operations (`src/fixed.rs`)

`saturating_add` / `saturating_sub` / `saturating_mul` return
`Saturated { value, saturated }` with the contract:

- `saturated == false` ⇒ `value` is the checked result: exact (add/sub)
  or within the multiplication floor bound (mul, Family C);
- `saturated == true` ⇒ `value` is clamped to `Fixed::MAX` / `Fixed::MIN`
  and **no numeric error bound is offered** — the exact result was out of
  range by an arbitrary amount. The flag is the signal to roll the phase
  back.

| Operation | Domain | Reverse operation | Class | Error |
|---|---|---|---|---|
| `saturating_add` | all `i64` inputs | — | lossy when saturated, exact otherwise | 0 if flag false; clamped (no bound) if true |
| `saturating_sub` | all `i64` inputs | — | same | same |
| `saturating_mul` | all `i64` inputs | — | same; unsaturated equals `checked_mul` | 0 within the Family C mul bound if flag false; clamped if true |

Tests: `tests/fixed_saturating_lossy.rs` — `saturating_add_sub_budget`
and `saturating_mul_budget_vs_i128_reference` (100k random pairs: the
flag is true iff the i128 reference is out of range, the value equals the
reference clamp, unsaturated equals the checked result).
`saturating_lossy_concrete_pairs` demonstrates the information loss:
`MAX.saturating_add(1)` and `MAX.saturating_add(2)` both return `MAX`
(the addend is unrecoverable), and the same is shown for sub and mul.

## Family E — `Interval<FRAC>` outward-rounded intervals (`src/interval.rs`)

Guarantee for every operation: for any real `x`, `y` inside the operand
intervals, the exact real result lies inside the returned interval — or
the operation returns `None` (an endpoint is unrepresentable, or the
divisor interval contains 0). Endpoints are rounded **outward**: `lo`
down (floor), `hi` up (ceiling). There is no infinity endpoint; refusal
(`None`) preserves soundness instead.

| Operation | Domain (`None` outside) | Class | Error / guarantee |
|---|---|---|---|
| `new(lo, hi)` | `lo ≤ hi` | exact | endpoints as given |
| `add(a, b)` | endpoints fit `i64` | exact, tight | containment exact: endpoints are the exact extremes, zero inflation |
| `sub(a, b)` | endpoints fit `i64` | exact, tight | containment exact: `[a.lo − b.hi, a.hi − b.lo]` |
| `mul(a, b)` | endpoints fit `i64` | per-endpoint lossy (outward) | contains every exact corner product; each endpoint is the outward rounding of the true extreme corner, hence within 1 ulp; width ≤ exact width + 2 ulp |
| `div(a, b)` | `b` excludes 0; endpoints fit `i64` | per-endpoint lossy (outward) | contains every exact corner quotient; same 1-ulp endpoint / +2-ulp width tightness |
| `sqrt(a)` | `lo ≥ 0`, `FRAC` even | per-endpoint lossy (outward) | contains every exact `sqrt(x)`; endpoints from floor/ceil sqrt, each within `2^−(FRAC/2)` |

Why corners suffice: `x·y` and `x/y` are monotone in each argument on
each side of 0, so over a box the extremes are attained at corners; the
divisor interval is required to exclude 0, so all four corner quotients
are finite and their min/max bound every quotient in the box. Floor is a
monotone map, so the minimum of the per-corner floors equals the floor of
the true minimum — the tests assert this as an equality, not just a bound.

Tests: `tests/interval_outward.rs` — `mul_containment_and_exact_outward_ends`
and `div_containment_and_exact_outward_ends` (50k random intervals each;
containment checked by sign-aware exact i128 cross-multiplication against
all four corners; endpoints asserted **equal** to the outward-rounded
exact extremes), `add_sub_exact_and_tight` (endpoints exact, width has
zero inflation, `None` ⇔ endpoint overflow), `sqrt_containment_and_outward_ends`,
`domain_errors_and_basics` (divisor containing 0, arithmetic overflow,
inverted bounds). The i128 floor/ceil division helpers are additionally
property-tested in `src/fixed.rs` (`floor_ceil_div_properties`).

## How the bounds are tested (claim → test map)

| Claim | Test |
|---|---|
| Family A inverses exact for every input | `tests/biject_inverse.rs` (all 5 tests) |
| Family C/D mul bound `−2^−F < err ≤ 0`, tight domain rule | `tests/fixed_bounds.rs::mul_floor_bound_f16`, `mul_floor_bound_f32` |
| Family C div bound (incl. negative divisors) | `tests/fixed_bounds.rs::div_floor_bound_f16` |
| Family C sqrt bound via true floor | `tests/fixed_bounds.rs::sqrt_bound_f16`, `sqrt_bound_f32` |
| Family C add/sub exact + inverse precondition | `tests/fixed_bounds.rs::add_sub_exact_and_invertible` |
| pow2/log2 exactness bands and roundtrip | `tests/fixed_bounds.rs::pow2_exactness_and_log2_roundtrip_f16`, `pow2_log2_f32`, `log2_floor_property` |
| Family B q_mul/q_div bounds + conservative domain | `tests/fixed_bounds.rs::legacy_q_mul_div_bounds_match_docs` |
| Family D saturating budget (flag ⇔ clamp) | `tests/fixed_saturating_lossy.rs::saturating_add_sub_budget`, `saturating_mul_budget_vs_i128_reference` |
| Lossy classes: concrete loss pairs | `tests/fixed_saturating_lossy.rs::saturating_lossy_concrete_pairs`, `mul_floor_lossy_concrete_pairs`, `sqrt_log2_lossy_concrete_pairs` |
| Partial inverses: documented preconditions | `tests/fixed_saturating_lossy.rs::partial_inverse_mul_div_precondition`, `partial_inverse_pow2_log2_precondition`, `fixed_add_sub_partial_inverse_documented`, `from_i32_inverse_precondition` |
| Family E containment + exact outward ends | `tests/interval_outward.rs::mul_containment_and_exact_outward_ends`, `div_containment_and_exact_outward_ends`, `add_sub_exact_and_tight`, `sqrt_containment_and_outward_ends` |
| Family E domain errors (no silent inflation) | `tests/interval_outward.rs::domain_errors_and_basics` |
| Exact decimal `Display` (no rounding) | `src/fixed.rs::tests::display_exact_decimal` |

All randomized tests use a deterministic xorshift64* generator with fixed
seeds, so the suite is reproducible in CI. Sample sizes: ≥ 50k per
operation family (200k for the mul/div bounds), ~1.7M operation checks in
total.

## Scope and honesty

- This document covers **implemented operations only**. There are no
  transcendentals beyond `2^x` (exact, or floor with error < 1 ulp below
  the resolution) and `floor(log2)`; no trigonometry, no exp/log beyond
  base 2, no floats, no NaN.
- The bounds are worst-case upper bounds over the stated domains,
  derived analytically and enforced by deterministic randomized tests
  against exact i128 references. They are test-verified guarantees, not
  machine-checked formal proofs.
- Intervals never inflate silently: if the true result cannot be
  represented, the operation returns `None`. Saturating arithmetic
  replaces a bound with a flag. In both cases the caller is expected to
  roll the phase back — that is the PHASE contract (fail = un-run).
- The legacy `q_*` functions are kept for compatibility and documented
  exactly as they behave, including the conservative `q_mul` domain and
  the truncating `q_div`.
