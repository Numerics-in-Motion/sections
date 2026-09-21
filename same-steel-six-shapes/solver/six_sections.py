# -*- coding: utf-8 -*-
"""Six sections holding the same steel -- the registered solver.

REGISTRATION frozen before any registered cell.
The existence gate ruled MAJOR REVISION on an earlier design; this file is
built to its five required changes and nothing here may be changed to suit a result.

THE QUESTION

An earlier study on this channel raced four sections holding the same 0.004000 m^2 of steel.
Two viewers asked where a triangular and a hexagonal hollow section would land. This adds
exactly those two.

WHAT IS HELD FIXED, AND WHY IT HAS TO BE STATED

Equal area does not fix a section's size -- that was that earlier study's whole point -- so a SECOND geometric
constraint is needed, and whichever one is chosen partly decides the answer. The registered rule
is **one common 200 x 200 mm bounding envelope**: every candidate is drawn as large as it fits
inside the same square hole, and its wall thickness is DERIVED from the area. One rule for the
whole roster. A per-shape outer ratio would be six hidden knobs, and that study's own config says
exactly that about its Rule B.

The rule is a CHOICE and is said on screen. An earlier design made the reversal between a round
and a square envelope its headline; the review rejected that, because an earlier study had
already established that changing the comparison rule can reverse a same-steel ranking.

PARENT AGREEMENT RUNS FIRST

`parent_agreement()` reproduces the earlier study's four published capacities before any new section is
computed, and `roster()` refuses to run if it fails. 014 shipped a result that reversed its own
parent's ranking because a section module was never called; that lesson is a gate here.
"""
from __future__ import annotations

import io
import json
import os

import numpy as np

import column_fem as CF
import column_geometry as CG
import section_properties as sp
from materials import Material

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG = os.path.join(ROOT, "config", "samesteel_config.json")

# --------------------------------------------------------------------------- registered constants
AREA = 0.004000              # m^2, the earlier study's steel area
ENVELOPE = 0.200             # m, the earlier study's envelope, used as a 200 x 200 mm bounding box
LENGTH = 3.000               # m, the earlier study's baseline -- NOT selected to clear a threshold
LENGTHS = (3.0, 4.0, 5.0, 6.0, 7.0, 8.0)     # the disclosed second axis
AREA_TOL = 1e-12             # m^2
PARENT_TOL = 1e-6            # relative, on the earlier study's four published capacities
TIE_FRACTION = 0.03          # the earlier study's own frozen tie fraction
I_FLANGE_RATIO = 0.5         # the earlier study's flange-width ratio
I_TF_TW = 2.0                # the earlier study's flange/web thickness ratio

# the earlier study's published capacities, quoted from its report -- never recomputed into this table
PARENT_PUBLISHED = {"solid_circle": 249233.02, "hollow_tube": 899508.44,
                    "box": 943972.83, "i_section": 356399.44}

NEW_POLYGONS = ((3, "triangle"), (6, "hexagon"))


def config():
    with io.open(CONFIG, encoding="utf-8") as f:
        return json.load(f)


def material():
    return Material.from_config(config()["material"])


# --------------------------------------------------------------------------- polygon geometry
def polygon_vertices(n, R):
    """A regular n-gon of circumradius R, FLAT SIDE DOWN.

    The orientation is registered, not incidental: a rotated polygon has a different bounding
    box and would be sized differently by the envelope rule.

    `pi/2 + pi/n` puts a flat side down only for EVEN n. For n = 3 it put a VERTEX down and a
    flat side up -- the opposite of the registration, and it shipped that way. An odd polygon
    needs no offset: a vertex at the top leaves a flat side at the bottom.
    """
    phase = np.pi / 2.0 + (np.pi / n if n % 2 == 0 else 0.0)
    th = phase + 2.0 * np.pi * np.arange(n) / n
    return np.c_[R * np.cos(th), R * np.sin(th)]


def shoelace(v):
    """Exact area and centroidal second moments of a simple polygon.

    Used instead of a quoted n-gon formula so that a mis-stated constant cannot enter: this
    integrates the actual vertices. `I_x == I_y` for a regular polygon, which `controls.py`
    checks rather than assumes.
    """
    x, y = v[:, 0], v[:, 1]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    cr = x * y1 - x1 * y
    a = cr.sum() / 2.0
    cx = ((x + x1) * cr).sum() / (6.0 * a)
    cy = ((y + y1) * cr).sum() / (6.0 * a)
    ixx = ((y * y + y * y1 + y1 * y1) * cr).sum() / 12.0 - a * cy * cy
    iyy = ((x * x + x * x1 + x1 * x1) * cr).sum() / 12.0 - a * cx * cx
    return abs(a), abs(ixx), abs(iyy)


def fit_circumradius(n, envelope=ENVELOPE):
    """The largest circumradius whose flat-side-down n-gon fits the square envelope."""
    v = polygon_vertices(n, 1.0)
    return envelope / max(np.ptp(v[:, 0]), np.ptp(v[:, 1]))


def hollow_polygon_properties(area, n, envelope=ENVELOPE):
    """Hollow regular n-gon carrying exactly `area`, sized to the envelope, wall DERIVED.

    The inner boundary is the outer polygon scaled about the centroid. Because the apothem
    scales with the circumradius, that is a constant wall thickness -- a real tube, not a
    mathematical convenience.
    """
    R = fit_circumradius(n, envelope)
    vo = polygon_vertices(n, R)
    a_out, ix_out, iy_out = shoelace(vo)
    if area >= a_out:
        raise ValueError("n=%d cannot carry %.6f m^2 inside a %.0f mm envelope"
                         % (n, area, 1000 * envelope))
    s = np.sqrt(1.0 - area / a_out)
    vi = polygon_vertices(n, R * s)
    a_in, ix_in, iy_in = shoelace(vi)
    apothem = R * np.cos(np.pi / n)
    ix, iy = ix_out - ix_in, iy_out - iy_in
    # The TRUE extreme fibre about each drawn axis. `c_x` goes with `I_x`, which is bending
    # about the horizontal axis, so it is the largest |y|. Setting both to R was geometrically
    # false: a flat-side-down hexagon reaches 0.0866 m vertically and 0.1000 m horizontally.
    # Equal inertia about every axis does not make the extreme fibre equal in every direction.
    cx_true = float(np.abs(vo[:, 1]).max())
    cy_true = float(np.abs(vo[:, 0]).max())
    return dict(kind="hollow_polygon_%d" % n, n=n, area=a_out - a_in, R=R, inner_scale=s,
                I_x=ix, I_y=iy, I_min=min(ix, iy), t=apothem * (1.0 - s),
                across_flats=2.0 * apothem,
                c_x=cx_true, c_y=cy_true,
                # the direction the COLUMN is analysed in -- see `solver_section`
                c_governing=R,
                outer_dim=2.0 * R,
                bbox_w=float(np.ptp(vo[:, 0])), bbox_h=float(np.ptp(vo[:, 1])),
                vertices_outer=vo, vertices_inner=vi)


# --------------------------------------------------------------------------- the roster
def sections(area=AREA, envelope=ENVELOPE):
    """The six registered sections, in the order they are drawn.

    the earlier four are rebuilt at the SAME envelope rule as the new two, so the comparison is one
    rule applied six times. Their published capacities are reproduced separately, by
    `parent_agreement()`, at that study's own proportions.
    """
    out = [("solid rod", sp.solid_circle_properties(area)),
           ("tube", sp.hollow_tube_properties(area, envelope)),
           ("box", sp.box_properties(area, envelope)),
           ("I-section", sp.i_section_properties(area, envelope,
                                                 I_FLANGE_RATIO * envelope, I_TF_TW))]
    for n, name in NEW_POLYGONS:
        out.append((name, hollow_polygon_properties(area, n, envelope)))
    for name, s in out:
        s.setdefault("I_min", min(s["I_x"], s["I_y"]))
        s["c_max"] = solver_section(s)["c_x"]        # the governing fibre, not the larger key
        err = abs(s["area"] - area)
        if err > AREA_TOL:
            raise ValueError("%s: area error %.3e exceeds %.1e" % (name, err, AREA_TOL))
    return out


def solver_section(sec):
    """The section AS THE SOLVER SEES IT: one inertia, and the extreme fibre that goes with it.

    A pin-ended column with no directional restraint may buckle in any direction. For the four
    parent sections the weak axis settles that, and `section_properties.fibre_for` returns the
    fibre belonging to `I_min`. For a regular polygon the inertia is the same about every
    centroidal axis, so the weak axis does NOT settle it -- what settles it is the extreme
    fibre, and the critical direction is the one that reaches furthest, which points at a
    VERTEX. That direction is registered and it is written down here rather than implied by a
    key someone has to guess at.
    """
    if "c_governing" in sec:
        c = sec["c_governing"]
    else:
        c = sp.fibre_for(sec, sec["I_min"])
    return dict(kind=sec["kind"], area=sec["area"],
                I_x=sec["I_min"], I_y=sec["I_min"], I_min=sec["I_min"],
                c_x=c, c_y=c, outer_dim=sec.get("outer_dim", 2.0 * c),
                c_governing=c)


# --------------------------------------------------------------------------- the parent gate
def parent_agreement(tol=PARENT_TOL):
    """Reproduce the earlier study's four published capacities. Runs BEFORE anything new is computed.

    an earlier study shipped a result that reversed its own parent's ranking because a section module was
    never called and every mass gate still passed. This is that lesson as a gate.
    """
    cfg = config()
    mat = material()
    s = cfg["solver"]
    rows = []
    for col in CG.build_all_columns(cfg):
        r = CF.run_nonlinear_column(col, mat, s["n_load_steps"], s["load_max_factor"],
                                    s["newton_max_iter"], s["newton_tol"])
        want = PARENT_PUBLISHED[col.name]
        got = float(r.failure_load)
        rows.append(dict(name=col.name, reproduced=got, published=want,
                         rel_err=abs(got - want) / want))
    worst = max(r["rel_err"] for r in rows)
    return dict(rows=rows, worst=worst, ok=bool(worst <= tol), tol=tol)


# --------------------------------------------------------------------------- capacity
def capacity(sec, L, mat=None):
    """THE REGISTERED ESTIMAND: first yield from the earlier study's own non-linear FEM.

    This called `analytical_imperfect_capacity` -- the closed form -- while the registration
    named `run_nonlinear_column`. The registration said one thing and the code did another, and
    nine gates passed because they all read this function's output. The difference is not
    cosmetic: the I-section moves about 15 %, and which pairs sit inside the 3 % tie rule
    changes with it.

    The closed form is still computed, by `capacity_closed_form`, as a CONTROL.
    """
    mat = mat or material()
    cfg = config()
    g, sv = cfg["geometry"], cfg["solver"]
    col = CG.build_column(sec["kind"], solver_section(sec), L, int(g["n_elements"]),
                          float(g["imperfection_ratio"]))
    r = CF.run_nonlinear_column(col, mat, sv["n_load_steps"], sv["load_max_factor"],
                                sv["newton_max_iter"], sv["newton_tol"])
    return float(r.failure_load)


def capacity_closed_form(sec, L, mat=None):
    """The parent's closed form, kept as an independent control on the FEM, never as the
    estimand."""
    mat = mat or material()
    ss = solver_section(sec)
    return float(CF.analytical_imperfect_capacity(
        mat.E, ss["I_min"], ss["area"], ss["c_x"], L,
        config()["geometry"]["imperfection_ratio"] * L, mat.failure_stress,
        config()["solver"]["effective_length_factor"]))


def euler(sec, L, mat=None):
    mat = mat or material()
    K = config()["solver"]["effective_length_factor"]
    return float(np.pi ** 2 * mat.E * sec["I_min"] / (K * L) ** 2)


def partial_order(rows, tie=TIE_FRACTION):
    """The order the result SUPPORTS: adjacent pairs inside the tie rule are not separated.

    A point-estimate rank is not a conclusion when the pair above is unresolved. This returns
    groups: consecutive sections whose gap is inside the rule share a group and are reported as
    unresolved with each other, never as first and second.
    """
    groups, cur = [], [rows[0]["name"]]
    for i in range(1, len(rows)):
        gap = 100.0 * (rows[i - 1]["P"] - rows[i]["P"]) / rows[i - 1]["P"]
        if gap < 100 * tie:
            cur.append(rows[i]["name"])
        else:
            groups.append(cur)
            cur = [rows[i]["name"]]
    groups.append(cur)
    return groups


def race(L=LENGTH, area=AREA, envelope=ENVELOPE):
    """The registered race at one length, ranked, with adjacent gaps."""
    mat = material()
    rows = []
    for name, s in sections(area, envelope):
        P = capacity(s, L, mat)
        rows.append(dict(name=name, P=P, Pcr=euler(s, L, mat), I=s["I_min"],
                         area=s["area"], c=s["c_max"], t=s.get("t"),
                         bbox_w=s.get("bbox_w", s.get("outer_dim")),
                         bbox_h=s.get("bbox_h", s.get("outer_dim")),
                         kind=s["kind"]))
    rows.sort(key=lambda r: -r["P"])
    for i, r in enumerate(rows):
        r["rank"] = i + 1
        r["gap_pct"] = (100.0 * (rows[i - 1]["P"] - r["P"]) / rows[i - 1]["P"]) if i else None
        r["tied_with_above"] = bool(i and r["gap_pct"] < 100 * TIE_FRACTION)
    return rows


def curve(lengths=LENGTHS, area=AREA, envelope=ENVELOPE):
    """The same six over the disclosed length axis -- shown, not selected from."""
    return [dict(L=L, rows=race(L, area, envelope)) for L in lengths]


__all__ = ["AREA", "ENVELOPE", "LENGTH", "LENGTHS", "TIE_FRACTION", "PARENT_PUBLISHED",
           "polygon_vertices", "shoelace", "fit_circumradius", "hollow_polygon_properties",
           "sections", "parent_agreement", "capacity", "euler", "race", "curve", "material"]
