# -*- coding: utf-8 -*-
"""The canonical run for this study. This is what ships.

    python src/six_sections_run.py

Writes `cells.json`: the parent-agreement gate, the registered race at the earlier study's own 3.000 m, the
disclosed length axis, and the kill line re-applied to the result.
"""
from __future__ import annotations

import io
import json
import os

import six_sections as S

DATA = os.environ.get("SIX_DATA") or os.path.abspath("six_data")


def build():
    parent = S.parent_agreement()
    if not parent["ok"]:
        raise SystemExit("PARENT AGREEMENT FAILED: worst %.3e > %.1e -- refusing to run"
                         % (parent["worst"], parent["tol"]))

    rows = S.race()
    by = {r["name"]: r for r in rows}
    tri, hexa = by["triangle"], by["hexagon"]
    tri_hex_gap = 100.0 * (hexa["P"] - tri["P"]) / hexa["P"]

    lengths = []
    for c in S.curve():
        b = {r["name"]: r for r in c["rows"]}
        lengths.append(dict(
            L=c["L"],
            order=[r["name"] for r in c["rows"]],
            P={r["name"]: r["P"] for r in c["rows"]},
            tri_hex_gap=100.0 * (b["hexagon"]["P"] - b["triangle"]["P"]) / b["hexagon"]["P"],
            spread=c["rows"][0]["P"] / c["rows"][-1]["P"]))

    fail = []
    if tri_hex_gap < 100 * S.TIE_FRACTION:
        fail.append("triangle-hexagon gap %.2f %% is inside the %.0f %% tie fraction"
                    % (tri_hex_gap, 100 * S.TIE_FRACTION))
    worst_area = max(abs(r["area"] - S.AREA) for r in rows)
    if worst_area > S.AREA_TOL:
        fail.append("area error %.3e exceeds %.1e" % (worst_area, S.AREA_TOL))

    res = dict(
        area=S.AREA, envelope=S.ENVELOPE, length=S.LENGTH,
        material=dict(E=S.material().E, fy=S.material().failure_stress,
                      rho=S.material().density),
        tie_fraction=S.TIE_FRACTION,
        parent=dict(worst=parent["worst"], tol=parent["tol"], rows=parent["rows"]),
        rows=[{k: v for k, v in r.items()} for r in rows],
        lengths=lengths,
        numbers=dict(
            tri_P=tri["P"], hex_P=hexa["P"], tri_rank=tri["rank"], hex_rank=hexa["rank"],
            tri_hex_gap=tri_hex_gap, hex_over_tri=hexa["P"] / tri["P"],
            spread=rows[0]["P"] / rows[-1]["P"],
            hex_over_i=hexa["P"] / by["I-section"]["P"],
            tri_over_i=tri["P"] / by["I-section"]["P"],
            hex_over_rod=hexa["P"] / by["solid rod"]["P"],
            tri_over_rod=tri["P"] / by["solid rod"]["P"],
            worst_area_err=worst_area),
        fail=fail, ship=not fail)

    os.makedirs(DATA, exist_ok=True)
    with io.open(os.path.join(DATA, "cells.json"), "w", encoding="utf-8", newline="\n") as f:
        json.dump(res, f, indent=1, sort_keys=True, default=float)
        f.write("\n")
    return res


def main():
    r = build()
    print("parent agreement worst %.3e  (tol %.0e)" % (r["parent"]["worst"], r["parent"]["tol"]))
    print("\nsame %.6f m2 of steel, %.0f x %.0f mm envelope, L = %.3f m\n"
          % (r["area"], 1000 * r["envelope"], 1000 * r["envelope"], r["length"]))
    for row in r["rows"]:
        g = "  --   " if row["gap_pct"] is None else "%5.2f %%" % row["gap_pct"]
        print("  %d. %-10s %8.1f kN   I %8.4f e-6 m4   gap %s%s"
              % (row["rank"], row["name"], row["P"] / 1e3, 1e6 * row["I"], g,
                 "  TIE" if row["tied_with_above"] else ""))
    n = r["numbers"]
    print("\n  hexagon rank %d, triangle rank %d, gap %.2f %% (kill line %.0f %%)"
          % (n["hex_rank"], n["tri_rank"], n["tri_hex_gap"], 100 * r["tie_fraction"]))
    print("  hexagon is %.2fx the I-section and %.2fx the rod; triangle %.2fx and %.2fx"
          % (n["hex_over_i"], n["hex_over_rod"], n["tri_over_i"], n["tri_over_rod"]))
    print("\n  length axis (disclosed, not selected from):")
    for L in r["lengths"]:
        print("    L = %.1f m  %-52s  gap %5.2f %%  spread %5.2fx"
              % (L["L"], " > ".join(L["order"]), L["tri_hex_gap"], L["spread"]))
    print("\n%s" % ("SHIPS" if r["ship"] else "KILLED: " + "; ".join(r["fail"])))
    print("-> %s" % os.path.join(DATA, "cells.json"))
    return 0 if r["ship"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
