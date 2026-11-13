"""Experimental Davies quadratic-form probabilities.

qf(chi2val, coeffs, dof=None, noncentrality=None, sigma=0.0,
   lim=1_000_000, acc=1e-8) returns (probability, fault, trace).

Numeric input sequences are extracted as f64; dof entries must be integers.
Faults follow Davies' diagnostics. Failed calculations return None for the
probability rather than the legacy wrapper's out-of-range 2.0 sentinel.
The legacy degenerate-boundary convention is retained: the probability at
statistic=0 for an identically zero form is one.
"""

from ._core import qf

__all__ = ["qf"]
