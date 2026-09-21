# Same steel, six shapes

Six cross-sections, each holding exactly **0.004000 m²** of steel, each drawn as large as it fits
inside the same **200 × 200 mm** envelope with its wall thickness derived. One pinned column,
**3.000 m**, with an L/1000 bow. The question two viewers asked of an earlier study: where do a
triangular and a hexagonal hollow section land?

**The triangle reaches first yield at 853.6 kN — fourth of six, and separated. The hexagon
reaches 911.0 kN, which is 6.30 % above the triangle but only 1.42 % below the tube: inside the
study’s pre-set 3 % pairwise tie rule, so it is not given a place of its own.**

The supported result is a PARTIAL order:

```
box  >  { tube , hexagon }  >  triangle  >  I-section  >  solid rod
             unresolved
```

In this model the first-yield loads of both new shapes are more than twice the I-section’s and
more than three times the solid rod’s.

*First-yield loads from an idealised planar-column model — not design capacities. The model
excludes local wall buckling and residual stresses; no design resistance is computed anywhere in
this repository.*

## The rule, which is a choice

Equal area does not fix a section's size — that was the earlier study's whole point — so a
second geometric constraint is needed, and whichever one is chosen partly decides the answer.
The rule here is **one common 200 × 200 mm bounding envelope**, applied to all six, with the
wall derived. One rule for the whole roster rather than a per-shape outer size, which would be
six hidden knobs. **A different rule can give a different order.**

Every polygon is drawn **flat side down**, and the extreme fibre is the vertex. Orientation is
part of the rule: a rotated polygon has a different bounding box and would be sized differently.

## Results at 3.000 m

`run_nonlinear_column`, the earlier study’s own non-linear FEM, is the estimand. The closed
form is a control and is not what is reported.

| rank | section | first yield (kN) | I_min (1e-6 m⁴) | gap to above (%) | |
|---|---|---|---|---|---|
| 1 | box | 958.4 | 25.3333 | - |  |
| 2 | tube | 924.1 | 18.7268 | 3.58 |  |
| 3 | hexagon | 911.0 | 15.3837 | 1.42 | tie |
| 4 | triangle | 853.6 | 11.7937 | 6.30 |  |
| 5 | I-section | 418.8 | 2.3358 | 50.94 |  |
| 6 | solid rod | 249.2 | 1.2732 | 40.49 |  |

A gap inside the pre-set 3 % rule means the pair is **not separated by this study’s
reporting convention** — not that the model failed to distinguish them. The rule is pairwise
and does not chain: box–tube is outside it, tube–hexagon inside it, box–hexagon outside it.

## The order does not depend on the length

| L (m) | order | triangle–hexagon gap (%) | spread |
|---|---|---|---|
| 3.0 | box > tube > hexagon > triangle > I-section > solid rod | 6.30 | 3.85 |
| 4.0 | box > tube > hexagon > triangle > I-section > solid rod | 10.45 | 6.36 |
| 5.0 | box > tube > hexagon > triangle > I-section > solid rod | 16.37 | 9.24 |
| 6.0 | box > tube > hexagon > triangle > I-section > solid rod | 20.63 | 12.30 |
| 7.0 | box > tube > hexagon > triangle > I-section > solid rod | 22.63 | 14.77 |
| 8.0 | box > tube > hexagon > triangle > I-section > solid rod | 23.34 | 16.35 |

The numerical order is identical at every length tested, and **every adjacent gap widens
monotonically** across the range. The length was not selected: 3.000 m is the earlier study’s
own baseline and is the first row. This is a robustness disclosure, not a second result.

## Parent agreement — the gate that runs first

Before any new section is computed, the solver reproduces the four capacities the earlier study
published, and **refuses to run if it cannot**:

| section | published (N) | reproduced (N) | rel. err |
|---|---|---|---|
| box | 943972.83 | 943972.83 | 1.82e-10 |
| hollow_tube | 899508.44 | 899508.44 | 2.36e-10 |
| i_section | 356399.44 | 356399.44 | 5.05e-09 |
| solid_circle | 249233.02 | 249233.02 | 1.13e-08 |

Worst 1.13e-08 against a tolerance of 1e-06. The four earlier shapes are reproduced **in their own
configuration**, not revised; the six shown above are recomputed under the common envelope.

## Kill line

Written before the registered run: the race dies if the triangle–hexagon gap — the comparison
that was actually asked for — falls inside the 3 % tie rule. It is **6.30 %**.

A clause was added after review: the race also dies if the triangle is not separated from
the group the hexagon sits in, because then there would be no placement to report at all.

## Reproduce

```
python solver/six_sections_run.py       # the parent gate, the race, the length axis
python solver/six_sections_checks.py    # C1-C8, the controls
```

`reference/reference_canonical.json` holds the canonical values. SHA-256 in that file is over
LF-normalised bytes.
