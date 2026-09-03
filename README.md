# Mesh to Part

A Fusion 360 add-in that turns an imported mesh (STL) into a **solid or a Form
body** under an explicit **maximum-error budget**, and then *measures* the real
deviation of the result against the original mesh instead of assuming it.

It includes a from-scratch organic mesh-to-Form conversion that does **not**
require the Product Design Extension.

![The add-in in use](docs/img/Mesh_control_panel.PNG)

| input mesh | Form body |
|---|---|
| ![Source mesh](docs/img/mesh-iso.png) | ![Form body](docs/img/form-iso.png) |

*A 10 876-triangle calcaneus scan, and the six-patch T-Spline Form body built from
it: watertight, volume within 0.23%, editable in the Form environment.*

## Why

Fusion's own **Convert Mesh → Organic** needs the Product Design Extension. Without
it you get `Prismatic`, which refuses organic anatomy outright, and `Faceted`,
which is exact but gives one BRep face per triangle. Neither is workable for a scan.

This add-in adds two methods of its own that need no extension:

| Method | Output | Needs extension |
|---|---|---|
| **Organic Form (T-Spline)** | a real **Form body**, 6 patches, editable with Fusion's sculpt tools | no |
| **Organic NURBS** | a 6-patch NURBS solid, watertight by construction | no |
| Prismatic | Fusion's own — planes, cylinders, cones | no |
| Faceted | Fusion's own — one face per triangle, exact | no |
| Organic | Fusion's own T-Spline conversion | **yes** |

Whatever the method, the budget is enforced the same way: convert, measure the
deviation **both ways** against the source mesh, and climb an escalation ladder
until it fits or the ladder runs out. The report says which happened.

## Measured accuracy

Every figure below was measured inside Fusion against the source mesh.

| input | method | max | RMS | p95 | volume error | faces |
|---|---|---|---|---|---|---|
| sphere ⌀30 mm | Form | 0.039 mm | 0.012 mm | 0.024 mm | +0.04% | 6 |
| sphere ⌀30 mm (engine) | NURBS | 0.012 mm | 0.004 mm | — | −0.02% | 6 |
| capsule (engine) | NURBS | 0.121 mm | 0.021 mm | — | −0.01% | 6 |
| calcaneus scan, 1063 mm diagonal | Form | 12.7 mm | 1.14 mm | 1.92 mm | −0.23% | 6 |
| calcaneus scan | NURBS | 12.9 mm | 0.87 mm | 1.28 mm | −0.13% | 6 |

A cage of radius 2.0 emitted as TSM produced a body **32.9146 mm** across; the
Catmull-Clark limit maths predicts **32.9146 mm**. A fitted 30 mm sphere came back
as **29.9986 mm**.

On the calcaneus, judge by p95/RMS. The max sits at a handful of isolated 150°
creases in a mesh whose mean edge length is 15.7 mm — no smooth surface reproduces
those, and Fusion's own Organic conversion would round them too.

## Install

1. Clone or copy this folder anywhere.
2. Run `powershell -File install.ps1`. It copies the add-in into
   `%APPDATA%\Autodesk\Autodesk Fusion 360\API\AddIns\MeshToPart`.
3. In Fusion: **Utilities → ADD-INS → Scripts and Add-Ins → Add-Ins** →
   `Mesh to Part` → tick **Run on Startup** → **Run**.
4. The button appears in **Mesh → Modify** and **Solid → Create**.

### The one external requirement

The organic fitting runs in a Python with **numpy, scipy and trimesh**, because
Fusion's bundled interpreter has none of them:

```
pip install numpy scipy trimesh
```

The add-in finds an interpreter automatically — the `MESHTOPART_PYTHON`
environment variable, a `python_path.txt` beside the add-in, the usual install
locations, then the `py` launcher. If none is found, the two organic methods are
simply absent from the method list. Prismatic and Faceted work without it.

## How the organic conversion works

A cube-sphere quad cage is shrink-wrapped onto the mesh, then relaxed
*tangentially* so points slide along the surface into concavities instead of being
lifted out of them. From there the two methods diverge:

- **Form** — the cage is solved so its Catmull-Clark *limit* surface lands on the
  mesh, then serialised as TSM and handed to `addByTSMDescription`.
- **NURBS** — each cube face becomes one bicubic patch fitted by constrained least
  squares; every shared cube edge is fitted **once** from the row of cage points
  both patches share, so the body is watertight by construction rather than by
  tolerance. `BRepBodyDefinition` assembles 8 vertices, 12 edges and 6 faces.

Full write-up, including what went wrong and what each mistake cost:
[docs/ORGANIC.md](docs/ORGANIC.md).

The TSM format is undocumented by Autodesk; what it is and how to write one is in
[docs/TSM-FORMAT.md](docs/TSM-FORMAT.md).

## Repository layout

```
MeshToPart.py     add-in shell: lifecycle, toolbar, events, dev bridge
mtp_core.py       all the working code; re-imported on every command
mtp_organic.py    the fitting engine, runs outside Fusion under numpy/scipy
install.ps1       push into Fusion's AddIns folder
resources/        toolbar icons
docs/             ORGANIC.md, TSM-FORMAT.md, index.html, screenshots
tools/            dev-bridge client and live test suites
tests/            offline tests, no Fusion required
```

`mtp_core.py` **hot-reloads** — it is re-imported on every command, so pushing it
takes effect on the next click; only `MeshToPart.py` needs a Stop/Run. Because of
that reload, nothing may keep state in a module global between events: the dialog's
method choice is resolved from the item label, never a remembered index.

## Tests

```
py -3.10 tests/test_geometry.py     # welding, gap spans, ladders, reporting
py -3.10 tests/test_engine.py       # cage, TSM structure, limit surface, fitting
```

`tests/stub` fakes `adsk`, so both run without Fusion. The engine tests check
against analytic answers — a 30 mm sphere, and the 32.9146 mm figure Fusion itself
produced from a known cage.

## The dev bridge

Fusion has no headless mode and no CLI, so the only way to test an add-in
automatically is from inside a running Fusion. Creating a `devbridge.on` file next
to the add-in starts a folder watcher that runs jobs on Fusion's main thread:

```
py -3.10 tools/mtp.py ping | probe | analyse | purge_work | rebuild_toolbar
py -3.10 tools/mtp.py convert '{"options": {"method_key": "form"}}'
py -3.10 tools/livetest.py      # full suite, everything in scratch documents
```

Off by default. Nothing in normal operation depends on it.

## Known limits

- The mesh-convert API is flagged *preview* by Autodesk; its enums can change.
- A Form body needs a **parametric** design. In a direct-modelling design the tool
  says so and points you at Organic NURBS instead.
- A finer cage follows a scan better but is likelier to self-intersect, which
  Fusion refuses outright. The Form ladder therefore steps *down* until a cage is
  accepted, and the report names the one it used.
- Preparation (weld, boundary walk) is pure Python inside Fusion — roughly 0.03 s
  per 10 k triangles. The organic fit is the slow part: cage 20 ≈ 30 s,
  cage 32 ≈ 70 s.

## Licence

MIT — see [LICENSE](LICENSE).
