# When plate buckling sets the first limit

One square hollow box. Hold the steel area at **0.004000 m²** and change nothing else: a thinner wall
buys a wider box, `B = A/4t + t`. Ask, at every thickness, which comes first -- elastic first
yield, or ideal plate buckling of the compression flange.

**For this elastic model, the moment at the first limit rises as the wall thins, peaks at
t = 4.6007 mm at 100.80 kN·m, and falls away below that.** The peak is where the two limits cross.

*This does not rate any real hollow section. Elastic first yield and an ideal plate limit only:
no residual stresses, no imperfections, no post-buckling reserve. No design resistance is
computed anywhere in this repository.*

## The model

| | |
|---|---|
| steel area | 0.004000 m², held exactly at every thickness |
| geometry | `B = A/(4t) + t`, flat width `b = A/(4t) - t` |
| first yield | `M_y = f_y Z_e`, f_y = 355 MPa (S355 to EN 10025) |
| ideal plate limit | `sigma_cr = k pi² E / (12 (1 - nu²)) (t/b)²`, k = 4, E = 210 GPa, nu = 0.30 |
| the estimand | `M_onset = Z_e min(f_y, sigma_cr)` |

`k = 4` is the classical coefficient for a long internal compression element with both edges
simply supported. It is the least favourable of the usual assumptions, and the sensitivity table
below carries it to 7.

## Results

| t (mm) | B (mm) | c/t | M_y (kN·m) | M_cr (kN·m) | governs |
|---|---|---|---|---|---|
| 12 | 95.3 | 5.94 | 35.2 | 2130.0 | yield |
| 10 | 110.0 | 9.00 | 43.5 | 1147.5 | yield |
| 8 | 133.0 | 14.62 | 55.8 | 558.3 | yield |
| 6 | 172.7 | 26.78 | 76.2 | 227.4 | yield |
| 5 | 205.0 | 39.00 | 92.4 | 129.9 | yield |
| 4 | 254.0 | 61.50 | 116.5 | 65.9 | ideal flange buckling |
| 3 | 336.3 | 110.11 | 156.4 | 27.6 | ideal flange buckling |

The seven rows are display cells. **The optimum is not one of them and is not found by searching
them:** it is the root

```
t_opt = sqrt(A / (4 (q + 1))),    q = sqrt(k pi² E / (12 (1 - nu²) f_y))
```

which gives **t = 4.6007 mm, B = 221.96 mm, c/t = 46.245, M_onset = 100.80 kN·m** -- between the 5 mm and
4 mm cells, and 9 % above the better of the two. A curve drawn by joining the seven cells
would put the maximum on the wrong one.

## Where this model leaves Class 3

`c/t = 42ε` (EN 1993-1-1 Table 5.2, internal compression part, ε = 0.8136) falls at
**t = 5.3321 mm**, where this model gives 86.39 kN·m. Below that thickness the section is Class 4.
The standard carries on from there through an effective cross-section; **this repository does
not, and computes no Class 4 design resistance anywhere.**

## Sensitivity to the plate coefficient

| k | t_opt (mm) | M_onset (kN·m) |
|---|---|---|
| 4 | 4.6007 | 100.80 |
| 5 | 4.3559 | 106.68 |
| 6 | 4.1653 | 111.73 |
| 7 | 4.0104 | 116.19 |

## Reproduce

```
python solver/box_wall_run.py        # the seven cells, the optimum, the Class 3 boundary
python solver/box_wall_checks.py     # C1-C9, the controls
```

`reference/reference_canonical.json` holds the canonical values. SHA-256 in that file is over
LF-normalised bytes.

## Why the controls ship with the solver

A closed-form model is easy to check and easy to believe without checking. Re-deriving the same
algebra a second time agrees to machine precision and establishes nothing, so each control is
answered a second way -- by integrating y² over the section rather than quoting `(B⁴ - Bi⁴)/12`,
by a golden-section search that knows nothing about the root, by the dimensional statement that
scaling A by λ² must scale every thickness by λ and every moment by λ³, and by the textbook the
plate constant came from.

Two of the nine were wrong when first written. C2 ran Simpson straight across the width
discontinuity at the inside face and was 1.1e-4 out -- indistinguishable from a wrong second
moment. C6's golden section mis-stated one of its four bracket updates and finished 47 % from a
root a plain grid finds to 5e-6 mm. Neither was caught by reading the code; both were caught by
running the controls, and each control is now also run against a model broken in exactly the way
it is meant to notice.
