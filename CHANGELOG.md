# Changelog

All notable changes to this project are documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versioning follows [SemVer](https://semver.org/).

## [Unreleased]

### Added

- Generic fixed-point type `Fixed<FRAC>` (value = `raw / 2^FRAC`, `i64` raw)
  with a documented Q-format scale: `checked_add`, `checked_sub`,
  `checked_neg`, `checked_mul`, `checked_div`, `sqrt` (even `FRAC`), `pow2`,
  `log2`, and an exact decimal `Display`; `Q16` alias for the classic
  Q16.16 layout. Every operation documents its domain and worst-case error
  bound; multiplication, division and sqrt round down (one-sided bounds).
- Saturating variants `saturating_add`, `saturating_sub`, `saturating_mul`
  returning `Saturated { value, saturated }`: the flag is `true` iff the
  result was clamped (no numeric bound when saturated; exact when not).
- Interval arithmetic `Interval<FRAC>` with outward rounding: `new`, `point`,
  `add`, `sub`, `mul`, `div`, `sqrt`, `contains`, `contains_zero`, `width`.
  Guarantee: the exact real result always stays inside the returned
  interval, or the operation returns `None` (no infinity endpoints);
  `mul`/`div` endpoints are exactly the outward roundings of the true
  extreme corners (within 1 ulp; width at most exact + 2 ulp).
- `docs/math-spec.md`: per-family specification tables (operation → domain →
  reverse operation → class → error bound / containment guarantee) with
  inline derivations and a claim-to-test map.
- Integration test suites (~1.7M operation checks, deterministic seeds):
  error bounds against exact i128 references (200k random full-range
  samples for the mul/div bounds, 100k for sqrt/log2), inverse-property
  tests for all bijective ops (edge values + random, all 64 rotations,
  ARX pairs, random pipelines), partial-inverse precondition tests,
  concrete lossy example pairs, saturating budget checks, interval
  containment/tightness tests.
- README quickstart (`cargo add` from git with a minimal example) and an
  honesty/scope section; Russian mirror `README.ru.md` with a language
  switcher.

### Changed

- Crate and API documentation translated to English; public API unchanged.
  The Q16.16 constants `FRAC` and `Q_MAX` are now re-exported at the crate
  root with documentation. Package version bumped to 0.2.0; crate metadata
  (repository, readme, keywords, categories) filled in.

### Fixed

- `q_pow2` documentation claimed exactness for the whole domain `|e| ≤ 31`;
  in fact it is exact only for `e ≥ −16` (for `e ≤ −17` it returns 0 with
  error `2^e < 1 ulp`). Behavior unchanged since v0.1.0; docs corrected and
  the exactness bands are now covered by tests.
- The shared Newton `isqrt` helper overflowed at `u64::MAX` (`x + 1` in the
  initial guess); it is now overflow-free over the whole `u64` range
  (previously latent and unreachable through the public API, which passes
  `raw ≤ i64::MAX`) and covered by unit tests.
