# Same steel, six shapes

Six cross-sections, each holding exactly **0.004000 m²** of steel, each drawn as large as it fits
inside the same **200 × 200 mm** envelope with its wall thickness derived. One pinned column,
**3.000 m**, with an L/1000 bow. The question two viewers asked of an earlier study: where do a
triangular and a hexagonal hollow section land?

**The hexagon lands third at 903.7 kN and the triangle fourth at 851.0 kN, 5.83 % apart.** Both
carry more than twice the I-section and more than three times the solid rod.

*The top three are within 3 % of the one above and this model does not separate them.*
*An idealised planar column at first yield: no local buckling of the walls, no residual
stresses, and no design resistance is computed anywhere in this repository.*

## The rule, which is a choice

Equal area does not fix a section's size — that was the earlier study's whole point — so a
second geometric constraint is needed, and whichever one is chosen partly decides the answer.
The rule here is **one common 200 × 200 mm bounding envelope**, applied to all six, with the
wall derived. One rule for the whole roster rather than a per-shape outer size, which would be
six hidden knobs. **A different rule can give a different order.**

Every polygon is drawn **flat side down**, and the extreme fibre is the vertex. Orientation is
part of the rule: a rotated polygon has a different bounding box and would be sized differently.

## Results at 3.000 m

| rank | section | first yield (kN) | I_min (1e-6 m⁴) | gap to above (%) | |
|---|---|---|---|---|---|
| 1 | box | 946.0 | 25.3333 | - |  |
| 2 | tube | 923.6 | 18.7268 | 2.36 | tie |
| 3 | hexagon | 903.7 | 15.3837 | 2.16 | tie |
| 4 | triangle | 851.0 | 11.7937 | 5.83 |  |
| 5 | I-section | 362.6 | 2.3358 | 57.39 |  |
| 6 | solid rod | 248.2 | 1.2732 | 31.54 |  |

## The order does not depend on the length

| L (m) | order | triangle–hexagon gap (%) | spread |
|---|---|---|---|
| 3.0 | box > tube > hexagon > triangle > I-section > solid rod | 5.83 | 3.81 |
| 4.0 | box > tube > hexagon > triangle > I-section > solid rod | 10.44 | 6.33 |
| 5.0 | box > tube > hexagon > triangle > I-section > solid rod | 16.24 | 9.26 |
| 6.0 | box > tube > hexagon > triangle > I-section > solid rod | 20.47 | 12.24 |
| 7.0 | box > tube > hexagon > triangle > I-section > solid rod | 22.47 | 14.72 |
| 8.0 | box > tube > hexagon > triangle > I-section > solid rod | 23.29 | 16.37 |

The ranking is identical at every length tested. Only the gaps widen. The length was not
selected: 3.000 m is the earlier study's own baseline and is the first row.

## Parent agreement — the gate that runs first

Before any new section is computed, the solver reproduces the four capacities the earlier study
published, and **refuses to run if it cannot**:

| section | published (N) | reproduced (N) | rel. err |
|---|---|---|---|
| box | 943972.83 | 943972.83 | 1.82e-10 |
| hollow_tube | 899508.44 | 899508.44 | 2.36e-10 |
| i_section | 356399.44 | 356399.44 | 5.05e-09 |
| solid_circle | 249233.02 | 249233.02 | 1.13e-08 |

Worst 1.13e-08 against a tolerance of 1e-06. The four earlier shapes are reproduced, not revised.

## Kill line

Written before the registered run: the race dies if the triangle–hexagon gap — the comparison
that was actually asked for — falls inside the 3 % tie fraction. It is **5.83 %**.

## Reproduce

```
python solver/six_sections_run.py       # the parent gate, the race, the length axis
python solver/six_sections_checks.py    # C1-C8, the controls
```

`reference/reference_canonical.json` holds the canonical values. SHA-256 in that file is over
LF-normalised bytes.
