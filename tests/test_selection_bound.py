import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from selection_bound import d_max_from_projections, d_max_many

rng = np.random.default_rng(0)


def test_no_loss_gives_zero():
    assert d_max_from_projections(rng.normal(size=1000), 1.0) == 0.0


def test_gaussian_matches_truncated_mean():
    # top S of a standard normal: E[Z | Z > z_{1-S}] = phi(z)/S
    from scipy.stats import norm
    p = rng.normal(size=400_000)
    for S in (0.5, 0.2, 0.05):
        z = norm.ppf(1 - S)
        assert abs(d_max_from_projections(p, S) - norm.pdf(z) / S) < 0.02


def test_bound_covers_any_selection_rule():
    # kill cells by an arbitrary state score; the induced mean shift along any u never exceeds d_max
    X = rng.normal(size=(20_000, 10)) * np.linspace(2, 0.5, 10)
    score = X[:, 0] + 0.5 * X[:, 3] + rng.normal(size=len(X))
    for S in (0.8, 0.4, 0.1):
        keep = score >= np.quantile(score, 1 - S)
        shift = X[keep].mean(0) - X.mean(0)
        u = shift / np.linalg.norm(shift)
        assert np.linalg.norm(shift) <= d_max_from_projections(X @ u, S) + 1e-9


def test_many_matches_single():
    X = rng.normal(size=(3000, 5)); U = rng.normal(size=(7, 5)); U /= np.linalg.norm(U, axis=1, keepdims=True)
    S = rng.uniform(0.05, 1, 7)
    many = d_max_many(X, U, S, chunk=3)
    single = [d_max_from_projections(X @ U[i], S[i]) for i in range(7)]
    assert np.allclose(many, single, atol=1e-4)
