# -*- coding: utf-8 -*-
"""Re-solve the registered box-wall cases and check every published number.

    python reproduce.py               # the seven cells, the optimum, the boundary + the controls
    python reproduce.py --quick       # the numbers only

Exits non-zero if any published value has moved, if the registered kill line fails, or if the
SHA-256 of a solver module no longer matches the frozen reference. Nothing here is fitted to the
answer: the reference was written from the canonical run before this file existed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "solver"))

REF = os.path.join(HERE, "reference", "reference_canonical.json")


def load_reference():
    with open(REF, encoding="utf-8") as f:
        return json.load(f)


def check_sha(ref):
    """SHA-256 over LF-normalised bytes, so a CRLF checkout still matches."""
    out = []
    for name, want in sorted(ref["source_sha256"].items()):
        p = os.path.join(HERE, "solver", name + ".py")
        with open(p, "rb") as f:
            got = hashlib.sha256(f.read().replace(b"\r\n", b"\n")).hexdigest()
        out.append(dict(module=name, ok=(got == want), got=got, want=want))
    return out


def solve_registered():
    import box_wall as BW
    cells = []
    for t_mm in BW.CELLS_MM:
        g = BW.onset(t_mm / 1000.0)
        cells.append(dict(t_mm=t_mm, B_mm=1000 * g["B"], c_over_t=g["c_over_t"],
                          My_kNm=g["My"] / 1e3, Mcr_kNm=g["Mcr"] / 1e3,
                          M_onset_kNm=g["M_onset"] / 1e3, governs=g["governs"]))
    t_opt = BW.t_optimum()
    go = BW.onset(t_opt)
    t_c3 = BW.t_at_slenderness(BW.class3_limit())
    gc = BW.onset(t_c3)
    return cells, dict(t_mm=1000 * t_opt, B_mm=1000 * go["B"], c_over_t=go["c_over_t"],
                       M_onset_kNm=go["M_onset"] / 1e3), \
        dict(t_mm=1000 * t_c3, M_onset_kNm=gc["M_onset"] / 1e3)


def compare(ref, cells, opt, c3):
    rel = ref["tolerances"]["moment_rel"]
    bad = []
    for row in cells:
        want = ref["primary"]["%.1f" % row["t_mm"]]
        for key in ("B_mm", "c_over_t", "My_kNm", "Mcr_kNm", "M_onset_kNm"):
            if abs(row[key] - want[key]) > rel * max(abs(want[key]), 1e-30):
                bad.append("t %.1f mm: %s %.10g, reference %.10g"
                           % (row["t_mm"], key, row[key], want[key]))
        if row["governs"] != want["governs"]:
            bad.append("t %.1f mm: governs %s, reference %s"
                       % (row["t_mm"], row["governs"], want["governs"]))
    for key, got, want in (("optimum t_mm", opt["t_mm"], ref["optimum"]["t_mm"]),
                           ("optimum M", opt["M_onset_kNm"], ref["optimum"]["M_onset_kNm"]),
                           ("class 3 t_mm", c3["t_mm"], ref["class3_boundary"]["t_mm"]),
                           ("class 3 M", c3["M_onset_kNm"],
                            ref["class3_boundary"]["M_onset_kNm"])):
        if abs(got - want) > rel * abs(want):
            bad.append("%s: %.10g, reference %.10g" % (key, got, want))
    return bad


def check_kill_line(ref, cells, opt):
    """The registered ship/kill line, re-applied rather than quoted.

    The study registered three things it could have failed on: the steel area must be held
    exactly at every thickness; the optimum must lie strictly INSIDE the displayed range, or
    there is no turn to show; and the limit that governs must actually change across it.
    """
    import box_wall as BW
    g = []
    worst = max(BW.geometry(c["t_mm"] / 1000.0)["area_err"] for c in cells)
    g.append(("steel area held to %.1e m^2" % worst, worst <= BW.AREA_TOL))
    ms = [c["M_onset_kNm"] for c in cells]
    g.append(("optimum inside the displayed cells",
              min(ms) < opt["M_onset_kNm"] > max(ms)))
    kinds = {c["governs"] for c in cells}
    g.append(("both limits govern somewhere", len(kinds) == 2))
    g.append(("the optimum is where they cross",
              abs(BW.onset(opt["t_mm"] / 1000.0)["My"]
                  - BW.onset(opt["t_mm"] / 1000.0)["Mcr"])
              <= 1e-9 * BW.onset(opt["t_mm"] / 1000.0)["My"]))
    return g


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="skip the controls")
    a = ap.parse_args(argv)
    ref = load_reference()
    ok = True

    print("SHA-256 of the solver modules (LF-normalised)")
    for r in check_sha(ref):
        print("  %-18s %s %s" % (r["module"], r["got"][:16], "OK" if r["ok"] else "CHANGED"))
        ok &= r["ok"]

    print("\nRe-solving the registered cases")
    cells, opt, c3 = solve_registered()
    bad = compare(ref, cells, opt, c3)
    for b in bad:
        print("  MOVED: %s" % b)
    ok &= not bad
    print("  %d cells, optimum t = %.4f mm at %.2f kN.m, class 3 boundary t = %.4f mm"
          % (len(cells), opt["t_mm"], opt["M_onset_kNm"], c3["t_mm"]))

    print("\nThe registered kill line")
    for name, good in check_kill_line(ref, cells, opt):
        print("  %-42s %s" % (name, "PASS" if good else "FAIL"))
        ok &= good

    if not a.quick:
        print("\nControls")
        r = subprocess.run([sys.executable, os.path.join(HERE, "solver", "box_wall_checks.py")],
                           capture_output=True, text=True)
        sys.stdout.write(r.stdout)
        ok &= (r.returncode == 0)

    print("\n%s" % ("EVERYTHING REPRODUCES" if ok else "SOMETHING HAS MOVED"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
