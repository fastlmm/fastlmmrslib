# Changelog

## 0.1.0 — Unreleased

- First Rust and Python release candidate, selected by the maintainer.
- Translate the full Davies quadratic-form calculation to safe Rust with
  per-call state; support signed weights, noncentrality, degrees of freedom,
  and an additional normal component.
- Return probability, Davies fault code, and seven diagnostic values. Failed
  probabilities are `None`, replacing the C++ wrapper's invalid 2.0 sentinel.
- Release the Python GIL during calculation. Test independent concurrent calls
  using Rayon in Rust and a thread pool through the Python binding.
- Detect non-finite parameter searches and cover the NaN-related stall regression.
- Qualify Python 3.11–3.15 installed wheels on Linux, macOS, and Windows, on
  x86-64 and ARM64, and source-derived wheels on boundary Python versions.

The prototype preceding this version passed all 38 native CI jobs. Final
0.1.0 artifacts require their own successful qualification before publication.
