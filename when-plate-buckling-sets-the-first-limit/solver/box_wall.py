# -*- coding: utf-8 -*-
"""How thin a box wall can get before local buckling arrives first -- the registered solver.

REGISTRATION (frozen by the design review before any registered cell was computed:
`analyze/260919/045_registration.md` v3, ruling MINOR REVISION). Nothing in this file may be
changed to suit a result.

THE QUESTION

An earlier study on this channel fixed the steel AREA at 0.004000 m^2 and showed that fixing
the area does not fix the SIZE. Thin
the wall at that fixed area and the box gets wider, so the elastic section modulus climbs -- until
the compression flange goes elastically critical before the steel yields. Where does that turn?

THE ESTIMAND, EXACTLY

    M_onset(t) = Z_e(t) * min( f_y , sigma_cr(t) )

the moment at the FIRST of elastic first yield or ideal elastic flange buckling. It is an ONSET,
not a capacity: a locally buckled plate keeps post-buckling reserve, and a compact section would
normally be assessed plastically. An earlier draft called it "the largest moment the section can
reach"; that is not supported and the review struck it.

WHAT THIS MODEL CANNOT SAY

Nothing about what a code would allow. EN 1993-1-1's class 4 is not where "the code stops you" --
it is where EN 1993-1-5's EFFECTIVE cross-section takes over, and no effective width is computed
here. A first draft claimed "the code gives up 15 %"; that was wrong and is deleted.

THE ANSWER IS SOLVED, NOT SEARCHED

sigma_cr = f_y  <=>  b/t = q = sqrt(k pi^2 E / (12 (1 - nu^2) f_y)), and with b = A/(4t) - t the
slenderness is b/t = A/(4 t^2) - 1, so

    t_opt = sqrt( A / (4 (q + 1)) )

exactly. The seven display thicknesses only draw the curve; they are not where the answer comes
from, and the curve itself is evaluated from the equations rather than interpolated between them.
"""
from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------- frozen registration
E = 210.0e9                  # Pa
NU = 0.30
FY = 355.0e6                 # Pa, S355 to EN 10025
AREA = 0.004000              # m^2 -- carried over unchanged from that earlier study
K_BASELINE = 4.0             # long internal compression element, both edges SIMPLY supported
K_SENSITIVITY = (4.0, 5.0, 6.0, 7.0)
CELLS_MM = (12.0, 10.0, 8.0, 6.0, 5.0, 4.0, 3.0)      # display cells: clean, round, 7 of them
AREA_TOL = 1e-12             # m^2, the geometry must keep the area exactly
T_MIN_MM, T_MAX_MM = 3.0, 12.0


def eps_steel(fy=FY):
    """EN 1993-1-1 epsilon = sqrt(235 / fy) with fy in MPa."""
    return float(np.sqrt(235.0e6 / fy))


def class3_limit(fy=FY):
    """Internal compression part, class 3 slenderness boundary: c/t = 42 eps."""
    return 42.0 * eps_steel(fy)


def geometry(t, area=AREA):
    """Square hollow box of wall t carrying exactly `area` of steel.

    A = B^2 - (B - 2t)^2 = 4 t B - 4 t^2  ->  B = A/(4t) + t, and the flange's clear width is
    b = B - 2t = A/(4t) - t.
    """
    t = float(t)
    B = area / (4.0 * t) + t
    b = area / (4.0 * t) - t
    if b <= 0.0:
        return None
    Bi = B - 2.0 * t
    I = (B ** 4 - Bi ** 4) / 12.0
    Ze = 2.0 * I / B
    return dict(t=t, B=B, Bi=Bi, b_flat=b, I=I, Ze=Ze, c_over_t=b / t,
                area_err=abs((B ** 2 - Bi ** 2) - area))


def sigma_cr(t, b, k=K_BASELINE):
    return k * np.pi ** 2 * E / (12.0 * (1.0 - NU ** 2)) * (t / b) ** 2


def onset(t, area=AREA, k=K_BASELINE):
    """The registered estimand at one thickness."""
    g = geometry(t, area)
    if g is None:
        return None
    scr = float(sigma_cr(g["t"], g["b_flat"], k))
    My = FY * g["Ze"]
    Mcr = scr * g["Ze"]
    g.update(k=k, sigma_cr=scr, My=float(My), Mcr=float(Mcr),
             M_onset=float(min(My, Mcr)),
             governs="yield" if My <= Mcr else "ideal flange buckling")
    return g


def q_ratio(k=K_BASELINE, fy=FY):
    """The slenderness at which sigma_cr = fy."""
    return float(np.sqrt(k * np.pi ** 2 * E / (12.0 * (1.0 - NU ** 2) * fy)))


def t_optimum(area=AREA, k=K_BASELINE, fy=FY):
    """The analytic root. No grid."""
    return float(np.sqrt(area / (4.0 * (q_ratio(k, fy) + 1.0))))


def t_at_slenderness(ratio, area=AREA):
    """t such that c/t is exactly `ratio` -- the same algebra, solved directly."""
    return float(np.sqrt(area / (4.0 * (float(ratio) + 1.0))))


def curve(n=600, t_min_mm=T_MIN_MM, t_max_mm=T_MAX_MM, area=AREA, k=K_BASELINE):
    """The curve for drawing, evaluated FROM THE EQUATIONS at every sample.

    The registered cells are seven points; joining them with straight lines badly understates how
    sharp the peak is, so the drawn curve is computed, not interpolated (the review's instruction).
    """
    ts = np.linspace(t_min_mm / 1000.0, t_max_mm / 1000.0, int(n))
    rows = [onset(t, area, k) for t in ts]
    return [r for r in rows if r]


__all__ = ["E", "NU", "FY", "AREA", "K_BASELINE", "K_SENSITIVITY", "CELLS_MM",
           "eps_steel", "class3_limit", "geometry", "sigma_cr", "onset",
           "q_ratio", "t_optimum", "t_at_slenderness", "curve"]
