# -*- coding: utf-8 -*-
"""Re-solve the six registered sections and check every published number.

    python reproduce.py               # the parent gate, the race, the length axis, the controls
    python reproduce.py --quick       # the numbers only

Exits non-zero if any published value has moved, if the registered kill line fails, or if the
SHA-256 of a solver module no longer matches the frozen reference.

The parent gate runs FIRST and the solver refuses to go on without it: before any new section is
computed, the four capacities an earlier study published are reproduced.
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


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="skip the controls")
    a = ap.parse_args(argv)
    ref = load_reference()
    ok = True

    print("SHA-256 of the solver modules (LF-normalised)")
    for r in check_sha(ref):
        print("  %-22s %s %s" % (r["module"], r["got"][:16], "OK" if r["ok"] else "CHANGED"))
        ok &= r["ok"]

    import six_sections as S

    print("\nParent agreement -- runs before any new section")
    pa = S.parent_agreement()
    for row in pa["rows"]:
        want = ref["parent_agreement"]["rows"][row["name"]]
        good = abs(row["reproduced"] - want["published_N"]) <= 1e-6 * want["published_N"]
        print("  %-14s %14.2f vs %14.2f  %.2e %s"
              % (row["name"], row["reproduced"], want["published_N"], row["rel_err"],
                 "OK" if good else "MOVED"))
        ok &= good
    ok &= pa["ok"]

    print("\nThe registered race at L = %.3f m" % ref["model"]["length_m"])
    rel = ref["tolerances"]["load_rel"]
    for row in S.race():
        want = ref["primary"][row["name"]]
        good = (row["rank"] == want["rank"]
                and abs(row["P"] - want["P_first_yield_N"]) <= rel * want["P_first_yield_N"])
        print("  %d. %-10s %9.1f kN  reference %9.1f kN  %s"
              % (row["rank"], row["name"], row["P"] / 1e3, want["P_first_yield_N"] / 1e3,
                 "OK" if good else "MOVED"))
        ok &= good

    print("\nThe kill line")
    rows = {r["name"]: r for r in S.race()}
    gap = 100.0 * (rows["hexagon"]["P"] - rows["triangle"]["P"]) / rows["hexagon"]["P"]
    good = gap >= 100 * ref["tie_fraction"]
    print("  triangle-hexagon gap %.2f %% against %.0f %%   %s"
          % (gap, 100 * ref["tie_fraction"], "PASS" if good else "FAIL"))
    ok &= good

    print("\nThe length axis -- one order at every length")
    orders = set()
    for L in ref["length_axis"]:
        got = [r["name"] for r in S.race(L["L_m"])]
        same = got == L["order"]
        orders.add(tuple(got))
        print("  L = %.1f m  %-52s %s" % (L["L_m"], " > ".join(got), "OK" if same else "MOVED"))
        ok &= same
    ok &= (len(orders) == 1)

    if not a.quick:
        print("\nControls")
        r = subprocess.run([sys.executable,
                            os.path.join(HERE, "solver", "six_sections_checks.py")],
                           capture_output=True, text=True,
                           env=dict(os.environ, PYTHONPATH=os.path.join(HERE, "solver")))
        sys.stdout.write(r.stdout)
        ok &= (r.returncode == 0)

    print("\n%s" % ("EVERYTHING REPRODUCES" if ok else "SOMETHING HAS MOVED"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
