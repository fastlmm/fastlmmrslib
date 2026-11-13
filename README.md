# fastlmmrslib

Experimental Rust successor to fastlmmclib's quadratic-form implementation.
This prototype is not a production replacement or a published release.
Version 0.0.0 is a development placeholder, not a maintainer-selected release.

The numerical core translates Robert Davies' full QFC algorithm into safe
Rust with per-call state and structured propagation of search-limit failures.
It supports signed coefficients, noncentrality, per-term degrees of freedom,
and an additional normal component. The Python binding releases the GIL.

```python
from fastlmmrslib import qf

probability, fault, trace = qf(20.0, [0.001, 0.001, 1.0, 1.0], acc=1e-7)
if fault != 0:
    raise ArithmeticError(f"Davies calculation reported fault {fault}: {trace}")
print(probability)  # approximately 4.5445e-5
```

The experimental Python interface retains the argument names used by
fastlmmclib. It returns `None` for a failed probability, rather than the legacy
wrapper's invalid value 2.0. Status 2 retains an approximate result with a
roundoff diagnostic; callers must inspect the status. Non-finite, mismatched,
or invalid numeric parameters return status 3. Type extraction errors raise
Python exceptions. This contract needs review before adoption.

The canonical Rust check is `cargo test`; strict linting is
`cargo clippy --all-targets -- -D warnings`. Build Python artifacts with
`maturin build --release --locked --interpreter PATH_TO_PYTHON`. Do not treat
a source-checkout test as installed-wheel qualification.

`compare_legacy.py` runs differential experiments against an installed
fastlmmclib artifact and the FaST-LMM evaluation inputs saved in
`/tmp/fastlmm-maintenance-review`. It is a development experiment, not a
standalone release gate; portable independent reference tests are still needed.

Retain `NOTICE.Davies` and `ThirdPartyNotices.original.txt` with distributions.
The original Microsoft-product/EULA preamble is retained for provenance,
not asserted as the licensing terms of this new project. Final provenance
wording requires human review before publication.
