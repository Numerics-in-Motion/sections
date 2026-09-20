# sections

Cross-sections of a fixed amount of material, solved from the equations and published with
everything needed to run them again.

| study | question |
|---|---|
| [`when-plate-buckling-sets-the-first-limit/`](when-plate-buckling-sets-the-first-limit/) | Hold the steel area and thin the wall of a square hollow box: the box gets wider and its elastic section modulus climbs, so how thin can the wall get before the compression flange goes elastically critical before the steel yields? |
| [`same-steel-six-shapes/`](same-steel-six-shapes/) | Six cross-sections holding the same steel area inside the same square envelope — a rod, a tube, a box, an I, a triangle and a hexagon. Where do the triangle and the hexagon land, and does the ranking depend on the column's length? |

Each study folder holds its solver, a frozen reference of every reported number with the SHA-256
of the modules that produced it, and a `reproduce.py` that re-solves the registered cases and
stops with an error if anything has moved.

```
pip install -r requirements.txt
cd when-plate-buckling-sets-the-first-limit
python reproduce.py --quick
```

Nothing here rates a real section. These are elastic models with their scope written down: first
yield and an ideal plate limit, no residual stresses, no imperfections, no post-buckling reserve,
and no design resistance is computed.

One of these studies reproduces the other's parent before it computes anything of its own: the
six-shape study re-runs an earlier study's four published capacities and **refuses to run if it
cannot**. A reader can check that from the folder.

The controls ship with the solver, and for a closed-form model they matter more rather than less:
re-deriving your own algebra a second time agrees to machine precision and establishes nothing.
Each control answers its question a second way, and two of the nine were wrong when first written
— one integrated straight across a discontinuity and was indistinguishable from a wrong second
moment, the other mis-stated a bracket update and finished 47 % from a root a plain grid finds to
five parts in a million. Neither was caught by reading the code.
