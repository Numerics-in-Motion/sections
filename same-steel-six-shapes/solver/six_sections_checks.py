# -*- coding: utf-8 -*-
"""Controls for the six-section study. C1-C8, run against the registered model.

    python src/six_sections_checks.py

Each answers its question a SECOND way. The polygon properties in particular are easy to get
wrong from a quoted formula, so nothing here quotes one: the areas and second moments come from
the shoelace integral over the actual vertices, and the controls test that integral against
cases whose answer is known independently.

C7 is the one that matters most for this subject. It tests the parent agreement itself against a
deliberately broken model, because a gate that has never rejected anything is a gate you are
trusting on its say-so.
"""
from __future__ import annotations

import numpy as np

import six_sections as S

TOL = dict(exact=1e-12, tight=1e-9)


def _report(name, ok, detail):
    print("  %-4s %-52s %s" % (name, detail, "PASS" if ok else "FAIL"))
    return ok


def c1_every_section_carries_the_registered_area():
    """The claim is "same steel". If the area drifted the comparison would be empty."""
    worst = max(abs(s["area"] - S.AREA) / S.AREA for _n, s in S.sections())
    return _report("C1", worst < TOL["exact"], "worst relative area error %.2e" % worst)


def c2_the_shoelace_integral_is_right_on_a_known_shape():
    """A square of side a has A = a^2 and I = a^4/12 about a centroidal axis. The polygon
    integral used for the triangle and the hexagon has to reproduce that exactly."""
    a = 0.137
    v = np.array([[-a / 2, -a / 2], [a / 2, -a / 2], [a / 2, a / 2], [-a / 2, a / 2]])
    A, ix, iy = S.shoelace(v)
    e = max(abs(A - a * a) / (a * a), abs(ix - a ** 4 / 12.0) / (a ** 4 / 12.0),
            abs(iy - a ** 4 / 12.0) / (a ** 4 / 12.0))
    return _report("C2", e < TOL["exact"], "square through the polygon integral %.2e" % e)


def c3_a_many_sided_polygon_converges_to_the_circle():
    """As n grows, a regular n-gon of circumradius R must approach pi R^2 and pi R^4 / 4.

    This is the only control that reaches outside the registered roster, and it is what says
    the integral is right for shapes it was not checked on."""
    R = 0.1
    A, ix, _iy = S.shoelace(S.polygon_vertices(4096, R))
    e = max(abs(A - np.pi * R ** 2) / (np.pi * R ** 2),
            abs(ix - np.pi * R ** 4 / 4.0) / (np.pi * R ** 4 / 4.0))
    return _report("C3", e < 1e-6, "n = 4096 against the circle %.2e" % e)


def c4_a_regular_polygon_has_no_weak_axis():
    """I is the same about every centroidal axis of a regular polygon, so the triangle and the
    hexagon have no orientation that is weaker in buckling. The renderer relies on this: it
    draws one orientation and reports `I_min`."""
    worst = 0.0
    for n in (3, 4, 5, 6, 8, 12):
        _a, ix, iy = S.shoelace(S.polygon_vertices(n, 0.1))
        worst = max(worst, abs(ix - iy) / ix)
        # and about a rotated axis, by rotating the shape instead
        th = 0.37
        rot = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
        _a2, ix2, _iy2 = S.shoelace(S.polygon_vertices(n, 0.1) @ rot.T)
        worst = max(worst, abs(ix2 - ix) / ix)
    return _report("C4", worst < 1e-9, "I about every axis, worst departure %.2e" % worst)


def c5_the_wall_is_a_constant_thickness():
    """The inner boundary is the outer polygon scaled about the centroid. That is a real tube
    only because the apothem scales with the circumradius; check the wall measured normal to a
    face equals the apothem difference."""
    worst = 0.0
    for n, _name in S.NEW_POLYGONS:
        s = S.hollow_polygon_properties(S.AREA, n)
        ap_out = s["R"] * np.cos(np.pi / n)
        ap_in = s["R"] * s["inner_scale"] * np.cos(np.pi / n)
        worst = max(worst, abs((ap_out - ap_in) - s["t"]) / s["t"])
    return _report("C5", worst < TOL["exact"], "wall vs apothem difference %.2e" % worst)


def c6_every_section_fits_the_envelope():
    """The rule is that each shape is as large as it FITS. Nothing may exceed the envelope."""
    worst = 0.0
    for name, s in S.sections():
        w = s.get("bbox_w", s.get("outer_dim"))
        h = s.get("bbox_h", s.get("outer_dim"))
        worst = max(worst, max(w, h) / S.ENVELOPE)
    return _report("C6", worst <= 1.0 + 1e-12, "largest bounding box / envelope %.6f" % worst)


def c7_the_parent_gate_rejects_a_broken_model():
    """The gate that reproduces the earlier study's four capacities is only worth something if
    it can fail. Break the material and require it to notice."""
    import column_fem as CF
    real = CF.run_nonlinear_column

    class _Bad:
        def __init__(self, r):
            self.failure_load = r.failure_load * 1.02
            self.__dict__.update({k: v for k, v in r.__dict__.items()
                                  if k != "failure_load"})

    CF.run_nonlinear_column = lambda *a, **k: _Bad(real(*a, **k))
    try:
        bad = S.parent_agreement()
    finally:
        CF.run_nonlinear_column = real
    good = S.parent_agreement()
    ok = (not bad["ok"]) and good["ok"]
    return _report("C7", ok, "2 %% error rejected (%.2e), truth accepted (%.2e)"
                   % (bad["worst"], good["worst"]))


def c8_the_capacity_is_bounded_by_squash_and_euler():
    """A first-yield load can exceed neither the squash load nor the Euler load. A closed form
    that violated either would be solving a different problem."""
    worst_s = worst_e = 0.0
    for name, s in S.sections():
        for L in S.LENGTHS:
            P = S.capacity(s, L)
            squash = s["area"] * S.material().failure_stress
            worst_s = max(worst_s, P / squash)
            worst_e = max(worst_e, P / S.euler(s, L))
    ok = worst_s <= 1.0 + 1e-9 and worst_e <= 1.0 + 1e-9
    return _report("C8", ok, "P/squash %.6f, P/Euler %.6f" % (worst_s, worst_e))


CHECKS = (c1_every_section_carries_the_registered_area,
          c2_the_shoelace_integral_is_right_on_a_known_shape,
          c3_a_many_sided_polygon_converges_to_the_circle,
          c4_a_regular_polygon_has_no_weak_axis,
          c5_the_wall_is_a_constant_thickness,
          c6_every_section_fits_the_envelope,
          c7_the_parent_gate_rejects_a_broken_model,
          c8_the_capacity_is_bounded_by_squash_and_euler)


def run():
    print("six sections controls")
    ok = [f() for f in CHECKS]
    print("  %d/%d pass" % (sum(ok), len(ok)))
    return all(ok)


if __name__ == "__main__":
    raise SystemExit(0 if run() else 1)
