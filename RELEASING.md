# Releasing fastlmmrslib

The maintainer selected 0.1.0 for the first release. Python and Rust versions
must agree in `pyproject.toml`, `Cargo.toml`, and `Cargo.lock`. No publication
or release tag is authorized merely by preparing this candidate.

## Qualify the candidate

- Run `cargo fmt --check`, `cargo test --locked`, and minimum-Rust 1.85 strict
  `cargo clippy --locked --all-targets --all-features -- -D warnings`.
  Also test on current stable Rust on every native platform. Build release
  wheels and sdist-derived wheels with current stable; 1.85 is a verified
  compatibility minimum, not the production compiler pin.
- Require the complete CI workflow: Rust checks on six native platforms,
  installed wheels on Python 3.11–3.15 on those platforms, and installed
  source-derived wheels on Python 3.11 and 3.15. The final job requires all
  30 wheels and one sdist and validates their metadata with pinned Twine.
  Linux release wheels are built and audited in manylinux 2.28 containers;
  require ten such wheels and reject host-specific `linux_*` wheel tags.
  Linux wheels therefore require glibc 2.28 or later. Local `uv build` wheels
  validate source/binding behavior but are not the Linux publication artifacts.
- Review independent analytic reference tests, the NaN-search stall regression,
  concurrent calls, the original notices, and `NOTICE.Davies` provenance.
  Retain original third-party notices verbatim. The historical Microsoft EULA
  preamble is preserved as provenance, not adopted as this project's license.
- Build wheel and sdist in isolation with `uv build --no-sources`; test the
  installed artifacts outside the checkout. Use the portable unittest suite
  in `tests/test_quadratic_form.py`. Differential comparison with the legacy
  C++ implementation is supplemental evidence, not the numerical oracle.
- Run `cargo publish --dry-run --locked` from a clean candidate commit.
  Check the crate's packaged sources, license, and notices.

## Maintainer publication

Before tagging, configure PyPI trusted publishing for owner `fastlmm`, repo
`fastlmmrslib`, workflow `release.yml`, and environment `pypi`; require maintainer
review on that environment. Review and approve the final provenance wording.
Merge the qualified candidate, date the changelog, and explicitly authorize
publication before creating or pushing `v0.1.0`.

The tag workflow verifies Python/Rust/tag version agreement and reuses the
complete CI workflow. Its publishing job waits for all qualification jobs and
the environment's approval, then publishes the exact
`fastlmmrslib-distributions` artifact. It never rebuilds after approval.
The Rust crate is separately published by the maintainer from the same clean
source with `cargo publish --locked`.

After publication, test clean installs of the published Python and Rust
packages. Only then change FaST-LMM's ordinary dependency graph and run final
acceptance against published fastlmmrslib 0.1.0. Do not qualify an ordinary
release by injecting a local sibling checkout or candidate wheel.

## Numerical integration contract

`qf` returns `(probability, fault, trace)`. Require fault zero before accepting
a FaST-LMM p-value. Fault 2 can include an approximate probability; other
failed calculations can return `None`. Do not turn these diagnostics into a
plausible p-value, clamp them, or fall back silently to the C++ backend.
The degenerate all-zero form retains the legacy boundary convention documented
in the Python module. Changing that convention requires a separate behavior
change and reference tests.
