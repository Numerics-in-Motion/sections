# -*- coding: utf-8 -*-
"""045 production run: the registered cells, the analytic event points, and the ship/kill line.

    python src/box_wall_run.py          # -> <data>/cells.json

Every gate in `analyze/260919/045_registration.md` §4 is re-applied here. The two answers -- the
optimum and the class-3 boundary -- come from the closed-form roots, never from the seven display
cells, so refining or coarsening the cells cannot move them.
"""
from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

import numpy as np                                            # noqa: E402

import box_wall as BW                                         # noqa: E402

# No internal numbering in anything that ships.
DATA = os.environ.get("BOX_DATA") or os.path.abspath("box_data")


def run():
    cells = []
    for t_mm in BW.CELLS_MM:
        r = BW.onset(t_mm / 1000.0)
        cells.append(dict(t_mm=t_mm, B_mm=1000 * r["B"], c_over_t=r["c_over_t"],
                          sigma_cr_MPa=r["sigma_cr"] / 1e6, My_kNm=r["My"] / 1e3,
                          Mcr_kNm=r["Mcr"] / 1e3, M_kNm=r["M_onset"] / 1e3,
                          governs=r["governs"], area_err=r["area_err"]))

    t_opt = BW.t_optimum()
    opt = BW.onset(t_opt)
    lim = BW.class3_limit()
    t_c3 = BW.t_at_slenderness(lim)
    c3 = BW.onset(t_c3)

    sens = []
    for k in BW.K_SENSITIVITY:
        tk = BW.t_optimum(k=k)
        rk = BW.onset(tk, k=k)
        sens.append(dict(k=k, t_mm=1000 * tk, B_mm=1000 * rk["B"], c_over_t=rk["c_over_t"],
                         M_kNm=rk["M_onset"] / 1e3))
    base = sens[0]["M_kNm"]
    for s in sens:
        s["pct"] = 100.0 * (s["M_kNm"] / base - 1.0)

    lo = [c for c in cells if c["t_mm"] == min(BW.CELLS_MM)][0]
    hi = [c for c in cells if c["t_mm"] == max(BW.CELLS_MM)][0]
    peak = opt["M_onset"] / 1e3

    fail = []
    if not (lo["M_kNm"] < peak > hi["M_kNm"]):
        fail.append("the optimum is not interior to the registered cells")
    if max(c["area_err"] for c in cells) > BW.AREA_TOL:
        fail.append("the geometry does not hold the registered steel area")
    if not (min(BW.CELLS_MM) < 1000 * t_opt < max(BW.CELLS_MM)):
        fail.append("the analytic optimum lies outside the displayed range")
    # every cell must be a valid section
    if any(c["c_over_t"] <= 0 for c in cells):
        fail.append("a cell has no flange left")

    out = dict(area=BW.AREA, fy=BW.FY, E=BW.E, nu=BW.NU, k=BW.K_BASELINE,
               eps=BW.eps_steel(), class3_limit=lim,
               cells=cells,
               optimum=dict(t_mm=1000 * t_opt, B_mm=1000 * opt["B"],
                            c_over_t=opt["c_over_t"], M_kNm=peak),
               class3=dict(t_mm=1000 * t_c3, c_over_t=c3["c_over_t"],
                           M_kNm=c3["M_onset"] / 1e3,
                           pct_of_peak=100.0 * c3["M_onset"] / opt["M_onset"]),
               ratios=dict(peak_over_thin=peak / lo["M_kNm"],
                           peak_over_thick=peak / hi["M_kNm"],
                           thin_pct=100.0 * lo["M_kNm"] / peak,
                           thick_pct=100.0 * hi["M_kNm"] / peak),
               sensitivity=sens, fail=fail, ship=not fail)
    os.makedirs(DATA, exist_ok=True)
    with open(os.path.join(DATA, "cells.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    return out


def main():
    o = run()
    print("Registered sweep: how thin a box wall can get before local buckling arrives first")
    print("square hollow box, steel area %.6f m2 fixed; S355 fy %.0f MPa, E %.0f GPa, k = %.1f"
          % (o["area"], o["fy"] / 1e6, o["E"] / 1e9, o["k"]))
    print()
    print("%-8s %-9s %-8s %-10s %-11s %-11s %-11s %s"
          % ("t (mm)", "B (mm)", "c/t", "sig_cr", "My (kNm)", "Mcr (kNm)", "M_onset", "governs"))
    for c in o["cells"]:
        print("%-8.1f %-9.1f %-8.1f %-10.0f %-11.2f %-11.2f %-11.2f %s"
              % (c["t_mm"], c["B_mm"], c["c_over_t"], c["sigma_cr_MPa"], c["My_kNm"],
                 c["Mcr_kNm"], c["M_kNm"], c["governs"]))
    print()
    op = o["optimum"]
    print("ANALYTIC optimum : t = %.4f mm, B = %.2f mm, c/t = %.3f, M_onset = %.2f kNm"
          % (op["t_mm"], op["B_mm"], op["c_over_t"], op["M_kNm"]))
    r = o["ratios"]
    print("   vs registered cells: %.1f mm is %.0f %% of the peak (%.2fx), %.1f mm is %.0f %% (%.2fx)"
          % (min(BW.CELLS_MM), r["thin_pct"], r["peak_over_thin"],
             max(BW.CELLS_MM), r["thick_pct"], r["peak_over_thick"]))
    c3 = o["class3"]
    print("EN class-3 bound : c/t = %.4f -> t = %.4f mm, M_onset = %.2f kNm (%.1f %% of the peak)"
          % (o["class3_limit"], c3["t_mm"], c3["M_kNm"], c3["pct_of_peak"]))
    print("   the ideal-plate optimum lies BEYOND that boundary, in class 4 territory;")
    print("   this model does not calculate class 4 design resistance.")
    print()
    print("sensitivity to the assumed edge restraint (baseline k = 4 is the reported case)")
    for s in o["sensitivity"]:
        print("   k %-4.1f  t %-7.2f mm  c/t %-6.1f  M_onset %-8.2f kNm  %+.1f %%"
              % (s["k"], s["t_mm"], s["c_over_t"], s["M_kNm"], s["pct"]))
    print()
    print("SHIPS" if o["ship"] else "KILLED: " + "; ".join(o["fail"]))
    return 0 if o["ship"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
