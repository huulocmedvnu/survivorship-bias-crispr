"""Trimming bound for survivorship (selection) bias in a perturbation direction.

If a treatment can only remove cells and the surviving fraction is S, the largest mean shift that selection alone
can produce along a unit direction u is

    d_max(S, u) = mean of the top ceil(S*n) control projections x.u  -  mean of all n control projections

(Lee 2009 trimming bound applied to the projection). Comparing d_max with the observed shift m along u gives an
upper bound on the share of the shift that selection can explain: ratio = d_max / m.
"""
import numpy as np


def d_max_from_projections(p, S):
    """p: 1-D array of control-cell projections on u (any centring); S: surviving fraction.
    Returns the largest mean shift along u that keeping a fraction S of the controls can produce."""
    p = np.asarray(p, dtype=np.float64)
    n = p.size
    if n == 0:
        return np.nan
    S = float(np.clip(S, 0.0, 1.0))
    k = max(1, int(np.ceil(S * n)))
    if k >= n:
        return 0.0
    top = np.partition(p, n - k)[n - k:]
    return float(top.mean() - p.mean())


def d_max_many(X, U, S, chunk=256):
    """X: (n_controls, d) control cells of one line; U: (k, d) unit directions; S: (k,) surviving fractions.
    Returns d_max for each direction. Works in chunks of directions to bound memory (n x chunk floats)."""
    X = np.asarray(X, dtype=np.float32)
    out = np.empty(len(U))
    for s in range(0, len(U), chunk):
        P = X @ np.asarray(U[s:s + chunk], dtype=np.float32).T          # (n, chunk)
        for j in range(P.shape[1]):
            out[s + j] = d_max_from_projections(P[:, j], S[s + j])
    return out
