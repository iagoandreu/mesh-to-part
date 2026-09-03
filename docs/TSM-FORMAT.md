# The TSM format

`TSplineBodies.addByTSMDescription()` and `addByTSMFile()` are the only routes in
the Fusion API to a T-Spline (Form) body, and the format they take is not
documented. The parser rejects everything malformed with a bare
`RuntimeError: invalid argument tsmDescription`, which gives you nothing to work
with.

This is what it actually is, reverse engineered from one hand-made Quad Ball
exported with `TSplineBody.getTSMDescription()`. A sample is in
[`quadball-sample.tsm`](quadball-sample.tsm); the writer is `write_tsm()` in
[`mtp_organic.py`](../mtp_organic.py).

It is **plain text**, not the binary container the `TSMZ` string in
`tsplines14.dll` suggests.

## Layout

```
#TS0200
                                  <- blank line
degree 3
cap-type G1CAPS
star-smoothness 0
units 1 meters
end-conditions SUBD_CREASES
star-knot-rule NURCCS
f <link> 0                        <- one per face
e <link> 1                        <- one per edge, 1 = smooth, 0 = creased
v <link> NORTH                    <- one per vertex
l <next> <prev> <mate> <vertex> <face> <edge> 0     <- one per half-edge
                                  <- blank line
0g <x> <y> <z> <w>                <- one per vertex, in v-record order
                                  <- blank line
tol 1.00000000000000008e-05
geom-tol 1.00000000000000008e-05
ver 1044
behavior-version 14.2.0
compat-version 13.3.0
```

It is a half-edge ("link") structure. `f`, `e` and `v` each name **one**
representative link; the rest of the connectivity is recovered by walking `l`.

## Verified invariants

Confirmed on Fusion's own Quad Ball (24 faces, 48 edges, 26 vertices, 96 links),
and asserted in [`tests/test_engine.py`](../tests/test_engine.py):

| field | meaning | check |
|---|---|---|
| `next` | next link around the face | `face[next[i]] == face[i]`, 4-cycles |
| `prev` | previous around the face | `prev[next[i]] == i` |
| `mate` | the same edge from the other face | `mate[mate[i]] == i`, `edge[mate[i]] == edge[i]`, `face[mate[i]] != face[i]` |
| `vertex` | the link's **destination** | `vertex[mate[i]] == vertex[prev[i]]` |
| `face` | face index | one `f` record per face |
| `edge` | edge index | exactly two links per edge |

**Winding.** The `next` loop runs **clockwise seen from outside** on all 24 faces
of the sample — the opposite of the usual counter-clockwise convention. A cage
built counter-clockwise has to be emitted reversed.

**Units.** Coordinates are in Fusion's internal unit, **centimetres**, despite the
`units 1 meters` line. Measured: a cage of radius 2.0 produced a body 32.9146 mm
across, and the Catmull-Clark limit radius of that cage is 1.6457 — so
1 unit = 10 mm exactly.

**`0g`** is `x y z w` homogeneous, `w = 1` for an unweighted control point. The
records are in the same order as the `v` records.

## Two things that will bite

**The body is the limit surface, not the cage.** A T-Spline body is the
Catmull-Clark limit surface of the control mesh you supply, and a subdivision
surface does not pass through its own control points. Emitting a cage that lies
on your target puts the body several millimetres inside it. The cage must be
solved for:

```
P_limit = (n² P + 4 Σ adjacent + Σ diagonal) / (n (n + 5))
```

so solve `L C = T` for the cage `C` given targets `T`. Relaxing toward it instead
diverges — the mask puts only `n/(n+5)` on the diagonal and the rest on the
neighbours, so Jacobi has no diagonal dominance to lean on. (Measured: the sphere
came out 536% over volume.)

**A self-intersecting cage produces no body at all.** Fusion rejects it with
`ASM_TSP_IG_SURFACE_PARALLEL_SINGULARITY` or `ASM_IG_SURFACE_SELF_INTERSECTS`, and
worse, the half-built form feature then poisons every later
`formFeatures.add()` with `InternalValidationError : baseFeature` until you delete
it. Solving `L C = T` exactly makes the control net overshoot in concavities — on
a real scan it turned 1–21 folded quads into 28–130 — so the solve needs a
regularisation term and the result needs a local untangle.

## Other API facts worth knowing

- `formFeatures.add()` and `FormFeature.startEdit()` work **without** the Product
  Design Extension. Only `Convert Mesh → Organic` is gated.
- The design must be **parametric**. In a direct-modelling design
  `formFeatures.add()` fails with `InternalValidationError : baseFeature`.
- Roll the timeline to the end (`design.timeline.moveToEnd()`) before adding a
  form feature. With the marker parked inside a base feature it fails with the
  same error.
- Fusion's own T-Spline primitive commands cannot be driven from the API to
  generate a sample: `commandCreated` does not fire for them to third-party
  handlers, and neither `ui.commandStarting` nor `ui.commandTerminated` hands over
  the `Command` object needed to press OK. One primitive made by hand is the only
  way in.
