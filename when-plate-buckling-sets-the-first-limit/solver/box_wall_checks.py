# -*- coding: utf-8 -*-
"""Controls for the box-wall study. C1-C9, run against the registered model.

    python src/box_wall_checks.py

WHY THESE SHIP WITH THE SOLVER

This model is closed form, which makes it easy to check and easy to believe without checking.
That is the trap: reproducing your own algebra, twice, agrees to machine precision and proves
nothing. So every control below is answered a SECOND way -- by numerical integration, by a
search that knows nothing about the root, by a dimensional argument, or by the textbook the
formula came from -- and the second way is never allowed to call the first.

The one that matters most is C2. `Ze` is `2I/B` with `I` written down in closed form; C2 rebuilds
`I` by integrating y^2 over the actual section, so a mis-stated second moment cannot survive.
"""
from __future__ import annotations

import numpy as np

import box_wall as BW

TOL = dict(exact=1e-12, tight=1e-9, loose=1e-6)
CELLS = [t / 1000.0 for t in BW.CELLS_MM]


def _report(name, ok, detail):
    print("  %-4s %-46s %s" % (name, detail, "PASS" if ok else "FAIL"))
    return ok


# ---------------------------------------------------------------- C1 geometry
def c1_area_is_conserved_exactly():
    """Every drawn section carries the registered steel area, to machine precision.

    The whole claim is "the same steel, every time". If the area drifted with t, the video would
    be comparing different amounts of steel and the comparison would be empty."""
    worst = 0.0
    for t in CELLS + [BW.t_optimum(), BW.t_at_slenderness(BW.class3_limit())]:
        g = BW.geometry(t)
        worst = max(worst, abs((g["B"] ** 2 - g["Bi"] ** 2) - BW.AREA) / BW.AREA)
    return _report("C1", worst < TOL["exact"], "area error %.2e (relative)" % worst)


# ---------------------------------------------------------------- C2 section properties
def c2_second_moment_by_integration():
    """Rebuild I by integrating y^2 dA over the real section, not from the closed form.

    The section is the region between two concentric squares, so at height y its width is B for
    |y| > Bi/2 and B - Bi inside. The only thing this control assumes is that shape; the closed
    form `(B^4 - Bi^4)/12` is never consulted, so a second moment written for the wrong shape
    cannot survive it.

    The width jumps at y = +-Bi/2, and the first version of this control ran Simpson straight
    across that jump on a uniform grid: 200 001 points and still 1.1e-4 out, which reads exactly
    like a wrong formula. The integral is split at the discontinuity instead.
    """
    from scipy.integrate import simpson
    worst = 0.0
    for t in CELLS:
        g = BW.geometry(t)
        B, Bi = g["B"], g["Bi"]
        I = 0.0
        for y0, y1, w in ((-B / 2.0, -Bi / 2.0, B), (-Bi / 2.0, Bi / 2.0, B - Bi),
                          (Bi / 2.0, B / 2.0, B)):
            y = np.linspace(y0, y1, 20001)
            I += w * float(simpson(y ** 2, x=y))
        worst = max(worst, abs(I - g["I"]) / g["I"])
    return _report("C2", worst < 1e-10, "I by integration differs by %.2e" % worst)


def c3_elastic_modulus_is_I_over_half_depth():
    worst = max(abs(BW.geometry(t)["Ze"] - 2.0 * BW.geometry(t)["I"] / BW.geometry(t)["B"])
                / BW.geometry(t)["Ze"] for t in CELLS)
    return _report("C3", worst < TOL["exact"], "Ze vs 2I/B differs by %.2e" % worst)


# ---------------------------------------------------------------- C4 the plate formula
def c4_critical_stress_scales_as_the_textbook_says():
    """sigma_cr must be exactly linear in k and exactly quadratic in t/b.

    The peak beat redraws the buckling limit for a moving k by ONE multiply on the strength of
    the first of those, so if it failed the video would be showing a curve nothing solved."""
    g = BW.geometry(0.005)
    s4 = BW.sigma_cr(g["t"], g["b_flat"], k=4.0)
    lin = max(abs(BW.sigma_cr(g["t"], g["b_flat"], k=k) - s4 * k / 4.0) / s4
              for k in (1.0, 4.0, 5.5, 7.0, 13.0))
    quad = abs(BW.sigma_cr(2 * g["t"], g["b_flat"]) - 4.0 * s4) / s4
    ok = lin < TOL["exact"] and quad < TOL["exact"]
    return _report("C4", ok, "k linearity %.2e, (t/b)^2 %.2e" % (lin, quad))


def c5_critical_stress_matches_the_classical_constant():
    """k pi^2 E / (12 (1 - nu^2)) written out independently, with the numbers substituted."""
    g = BW.geometry(0.006)
    want = (4.0 * (np.pi ** 2) * 210.0e9 / (12.0 * (1.0 - 0.30 ** 2))
            * (g["t"] / g["b_flat"]) ** 2)
    rel = abs(BW.sigma_cr(g["t"], g["b_flat"]) - want) / want
    return _report("C5", rel < TOL["exact"], "sigma_cr vs classical form %.2e" % rel)


# ---------------------------------------------------------------- C6 the optimum
def c6_the_optimum_is_found_by_a_search_that_does_not_know_the_root():
    """Golden-section on M_onset, started from the ends of the drawn range.

    `t_optimum` is an analytic root. A search agreeing with it is the only evidence that the
    root is the root of the right function."""
    a, b = 0.0025, 0.0110
    phi = (np.sqrt(5.0) - 1.0) / 2.0
    f = lambda t: BW.onset(t)["M_onset"]
    # the bracket is recomputed each pass rather than carried forward: the first version tried
    # to reuse the interior points and mis-stated one of the four updates, which moved the
    # bracket the wrong way and landed 47 % off a root a plain grid finds to 5e-6 mm
    for _ in range(120):
        c, d = b - phi * (b - a), a + phi * (b - a)
        if f(c) > f(d):
            b = d
        else:
            a = c
    found = 0.5 * (a + b)
    rel = abs(found - BW.t_optimum()) / BW.t_optimum()
    return _report("C6", rel < 1e-7, "golden section vs analytic root %.2e" % rel)


def c7_at_the_optimum_the_two_limits_are_equal():
    """The peak is a CROSSING, not a smooth maximum: it is where first yield and ideal plate
    buckling coincide. If they did not, the closed form would be solving a different problem."""
    g = BW.onset(BW.t_optimum())
    rel = abs(g["My"] - g["Mcr"]) / g["My"]
    return _report("C7", rel < 1e-10, "M_y vs M_cr at the optimum %.2e" % rel)


# ---------------------------------------------------------------- C8 the code boundary
def c8_the_class_3_boundary_is_42_epsilon():
    """EN 1993-1-1 Table 5.2, internal compression part: c/t = 42 eps, eps = sqrt(235/f_y)."""
    eps = np.sqrt(235.0 / 355.0)
    lim = 42.0 * eps
    t = BW.t_at_slenderness(lim)
    g = BW.geometry(t)
    ok = (abs(BW.class3_limit() - lim) < TOL["exact"]
          and abs(g["c_over_t"] - lim) / lim < TOL["tight"])
    return _report("C8", ok, "c/t at the boundary %.5f vs 42 eps %.5f" % (g["c_over_t"], lim))


# ---------------------------------------------------------------- C9 similarity
def c9_the_model_scales_like_a_length():
    """Scale the area by lambda^2 with f_y fixed: every thickness must scale by lambda and every
    moment by lambda^3. A control the closed form cannot pass by accident, because it mixes the
    geometry, the stress limit and the optimum in one statement."""
    worst = 0.0
    for lam in (0.5, 1.37, 2.0, 3.3):
        A2 = BW.AREA * lam ** 2
        worst = max(worst, abs(BW.t_optimum(area=A2) - lam * BW.t_optimum()) / BW.t_optimum())
        m1 = BW.onset(BW.t_optimum(), BW.AREA)["M_onset"]
        m2 = BW.onset(BW.t_optimum(area=A2), A2)["M_onset"]
        worst = max(worst, abs(m2 - (lam ** 3) * m1) / ((lam ** 3) * m1))
    return _report("C9", worst < TOL["tight"], "similarity residual %.2e" % worst)


CHECKS = (c1_area_is_conserved_exactly, c2_second_moment_by_integration,
          c3_elastic_modulus_is_I_over_half_depth,
          c4_critical_stress_scales_as_the_textbook_says,
          c5_critical_stress_matches_the_classical_constant,
          c6_the_optimum_is_found_by_a_search_that_does_not_know_the_root,
          c7_at_the_optimum_the_two_limits_are_equal,
          c8_the_class_3_boundary_is_42_epsilon,
          c9_the_model_scales_like_a_length)


def run():
    print("box wall controls")
    ok = [f() for f in CHECKS]
    print("  %d/%d pass" % (sum(ok), len(ok)))
    return all(ok)


if __name__ == "__main__":
    raise SystemExit(0 if run() else 1)
