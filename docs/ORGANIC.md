# The organic conversion

Fusion's **Convert Mesh → Organic** requires the Product Design Extension. This is
a from-scratch replacement that produces either a real T-Spline **Form body** or a
smooth six-patch **NURBS solid** from a mesh, using only APIs the base product
exposes.

## Two outputs

**Organic Form (T-Spline)** is the real thing: a Form body you can sculpt with
Fusion's own tools. It works because `formFeatures.add()` and `startEdit()` are
*not* gated by the Product Design Extension — only `Convert Mesh → Organic` is.
The missing piece was the TSM format, which Autodesk does not document and whose
parser answers every mistake with a bare `invalid argument tsmDescription`. It was
reverse engineered from a single hand-made Quad Ball; see
[TSM-FORMAT.md](TSM-FORMAT.md).

**Organic NURBS** builds a six-patch NURBS solid through `BRepBodyDefinition`
instead. It is more accurate on a difficult scan, works in direct-modelling
designs where Form features are unavailable, and does not care about
self-intersection the way the T-Spline converter does.

Both share the cage machinery below; they differ only in what is done with the
fitted cage.

## The method

**1. Cube-sphere cage.** Six `n × n` quad grids on a cube, projected to a sphere, with
shared cube edges and corners resolving to the same vertex. No pole singularity and no
periodic seam, so every patch is a plain open bicubic surface and the topology is the
simplest watertight quad layout that exists: 6 faces, 12 edges, 8 vertices.

**2. Shrink-wrap.** The cage starts as a sphere enclosing the mesh and is pulled in
against a KD-tree over a dense surface sample, with Laplacian smoothing to stop it
folding on the way.

**3. Tangential relaxation.** Every cage point is projected exactly onto the mesh and
then only allowed to slide *along* it — the Laplacian displacement is stripped of its
normal component before it is applied. This is the step that reaches into concavities:
ordinary smoothing lifts a point out of a hollow and the next projection snaps it back
to the nearest convex face, so the hollow never fills. Normals are averaged over the
cage neighbourhood first, because a raw face normal flips across a crease and makes the
relaxation oscillate instead of settle; steps are capped at 0.35 of local spacing.

Measured result: the cage lands on the mesh to **1e-12 mm**.

**4. Constrained least-squares patch fitting.** Each cube face becomes one bicubic
NURBS patch. Every cube edge is fitted **once**, from the row of cage points the two
adjacent patches share, with both end control points pinned to the corner points. Both
patches then use that identical boundary curve, so the body is watertight *by
construction* rather than by tolerance — and three patches meeting at a corner agree on
that corner exactly. The interior control points solve a least-squares problem with the
boundary rows held fixed and a second-difference fairness term.

Two things this got wrong first, both caught by measurement:
- Pure interpolation is watertight and passes through every cage point, but rings
  between them — 115 of 98 304 sampled points were more than 5 mm out, in patch
  interiors, never at boundaries.
- Solving via the normal equations (`AᵀA`) squares the condition number, which at full
  control density put **114 mm** of ringing into the fit. The stacked system is solved
  directly with `lstsq` instead.

**5. B-Rep assembly.** `BRepBodyDefinition`: 8 vertex definitions, 12 edge definitions
carrying the shared `NurbsCurve3D`, 6 face definitions carrying the `NurbsSurface`, each
with one loop of 4 co-edges. `CUBE_FACES` is built so `u × v` is the outward normal,
which makes the parameter-space loop counter-clockwise from outside and means no face
needs `isParamReversed`.

Two API facts established by test, not by guessing: `NurbsSurface.create` wants an
**empty** weights array for a non-rational surface (a list of 1.0s is rejected with
`invalid argument weights`), and the flat control point array is **u-major** — probed by
building a deliberately asymmetric 4×5 net and evaluating a known corner.

## Where the work runs

Fusion's bundled interpreter has no numpy, so the fitting runs in `py -3.10`
(numpy 2.2.6 / scipy 1.15.3 / trimesh 4.12.2), launched as a subprocess from the
add-in. `rtree` is not installed and was not added: exact closest-point queries use a
KD-tree over triangle centroids to shortlist candidates, bounded so the answer is exact,
plus a vectorised point-triangle projection.

## Measured accuracy

Deviation is always reported both ways and judged on the worse one.

| input | cage | part → mesh max | RMS | p95 | mesh → part max | folds | volume err |
|---|---|---|---|---|---|---|---|
| icosphere ⌀30 mm | 16 | **0.012 mm** | **0.004 mm** | — | 0.017 mm | 0 | −0.024% |
| capsule | 20 | 0.121 mm | 0.021 mm | — | 0.049 mm | 0 | −0.008% |
| sphere ⌀30 mm, in Fusion | 16 | 0.039 mm | 0.012 mm | 0.024 mm | 0.049 mm | 0 | +0.04% |
| calcaneus scan, 1063 mm diag | 32 | 12.9 mm | 0.87 mm | **1.28 mm** | 9.4 mm | 64 | −0.128% |
| calcaneus scan | 48 | 17.3 mm | 0.70 mm | **0.88 mm** | 6.4 mm | — | −0.16% |

All rows are watertight 6-face solids. The calcaneus rows are measured by Fusion
against the original mesh, independently of the engine's own check.

On clean geometry the method is essentially exact. The calcaneus is the hard case, and
its numbers are dominated by the input rather than the fit: mean edge length 15.7 mm
(max 354 mm), normals disagreeing by up to 153° across creases, and concavities 110 mm
deep. RMS, p95 and p99 improve steadily with cage resolution (engine RMS 0.60 → 0.48 →
0.43 mm at cage 32 → 40 → 48; p95 1.28 → 0.88 mm) while the max stays around 13-17 mm
at a handful of isolated spots — a smooth surface cannot reproduce a 150° crease, and
neither could Fusion's own Organic conversion. Judge this method on p95/RMS; the max
is a crease report, not a fit quality report.

Note that mesh is at roughly 10× scale for a calcaneus; at true size those figures are
RMS 0.04 mm and max 1.3 mm.

## The Form path

The cage is shared with the NURBS path up to the point where it lies on the mesh.
Then:

1. **Solve for the limit surface.** A T-Spline body is the Catmull-Clark limit
   surface of its control cage, and that surface does not pass through the control
   points — a cage lying on the scan gives a body several millimetres inside it. So
   the cage is solved for, `L C = T`. Relaxing toward it instead diverges: the
   limit mask puts only `n/(n+5)` on the diagonal, leaving Jacobi no diagonal
   dominance, and the sphere came out **536% over volume**.
2. **Regularise, then untangle locally.** Solving `L C = T` exactly makes the
   control net overshoot in concavities: on the calcaneus it turned 1–21 folded
   quads into 28–130. Raising a global smoothing term until folds cleared worked
   but cost a factor of four in RMS. A light term plus a *local* untangle of the
   folded neighbourhoods keeps the accuracy the fit earned.
3. **Emit TSM and let Fusion build it.** Folds are fatal here, not cosmetic:
   Fusion refuses a self-intersecting T-Spline outright and returns no body. A
   finer cage follows the scan better but is likelier to self-intersect, so the
   Form ladder steps **down** until a cage is accepted — on the calcaneus, cage 28
   was refused and cage 21 accepted.

A rejected cage also leaves the form feature half built, after which every later
`formFeatures.add()` fails with `InternalValidationError : baseFeature`; tearing it
down on failure is what lets the ladder try the next rung at all.

## Cage folding — the remaining limitation

On a strongly concave shape the cage can buckle: a quad ends up facing away from the
mesh underneath it, and the fitted surface self-intersects there. Three attempts, in
order, and what each taught:

1. **Repair afterwards** (full Laplacian on folded regions, reproject). Stalled at
   65 of 75 folds and made the fit *worse* — RMS 0.73 → 0.81 mm.
2. **Prevent per step, self-consistency test** (did this quad's normal flip during
   this step?). Caught nothing: the count climbed steadily from 1 to 80 during the
   wrap while no single step ever flipped a quad. Buckles build up gradually.
3. **Prevent per step, against the mesh normal.** Correct test — but the first
   version also backed off around quads that were *already* folded, which starved
   those vertices of movement and left the cage 72 mm short of the surface in
   places. Restricting the back-off to folds the step itself creates fixed that.

Result: folds after the wrap fell from 139 to 85 on the calcaneus, and the wrap
residual returned to 2.7 mm. About 2% of quads still fold on that mesh and resist
further relaxation. **On clean geometry there are none** — the sphere finishes with
zero folded quads and 0% inverted facets.

The residue is where a genus-0 cube-sphere parameterisation is genuinely strained: a
single quad chart has to cover a bone with 110 mm-deep concavities. Fixing it properly
means a feature-aligned or multi-chart layout, which is a different algorithm rather
than a tuning change. The fold count is reported in every run so it is never silent.

## Controls

- **Maximum error** drives everything. With *cage from the error budget* ticked, the
  starting cage is `0.06 × bbox_diagonal / max_error`, clamped to 20…48, and the
  ladder then climbs `cage × 1.35` per attempt against the measurement. The floor is
  20 because below that the surface is too coarse to be worth building whatever the
  budget says — a first version with a floor of 12 dutifully produced a cage-12 body
  at 2.3 mm RMS when asked for 10 mm.
- **Fairness** trades ringing against fidelity; 0.01 is the tested default.
- **Relax passes** — more is better on concave shapes and costs roughly linear time.

Cost scales steeply: cage 20 ≈ 30 s, cage 32 ≈ 70 s, cage 48 ≈ 160 s for the engine
alone, plus the Fusion-side deviation measurement.
