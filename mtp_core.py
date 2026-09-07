"""Mesh to Part - all the working code. Re-imported fresh on every use.

MeshToPart.py is the shell (lifecycle, toolbar, event wiring). Everything that is
worth iterating on lives here, so editing this file takes effect on the next command
or the next dev job - no Stop/Run, no Fusion restart.

Internal Fusion length unit is CENTIMETRES throughout. Everything the user sees is
millimetres; the conversion happens only at the input and report boundaries.
Fusion's bundled Python has no numpy - stdlib only in this file.
"""

import adsk.core
import adsk.fusion

import json
import math
import os
import sys
import time
import traceback
import types

VERSION = '2.0.0'

ADDIN_DIR = os.path.dirname(os.path.realpath(__file__))
SETTINGS_PATH = os.path.join(ADDIN_DIR, 'last_settings.json')
LOG_PATH = os.path.join(ADDIN_DIR, 'MeshToPart.log')
REPORT_DIR = os.path.join(ADDIN_DIR, 'reports')

CM = 10.0          # cm -> mm
MM = 0.1           # mm -> cm
WELD_TOL_CM = 1.0e-7


def log(msg):
    try:
        with open(LOG_PATH, 'a', encoding='utf-8') as fh:
            fh.write('%s  %s\n' % (time.strftime('%H:%M:%S'), msg))
    except Exception:
        pass


def short_error(exc_text):
    """Last meaningful line of a traceback, for a human-facing message."""
    lines = [l.strip() for l in exc_text.strip().splitlines() if l.strip()]
    for l in reversed(lines):
        if not l.startswith('File "') and not l.startswith('~') and not l.startswith('^'):
            return l
    return lines[-1] if lines else 'unknown error'


# --------------------------------------------------------------------------- #
# enumerations
# --------------------------------------------------------------------------- #

ACC_TIERS = [
    ('Low', adsk.fusion.MeshConvertAccuracyTypes.LowMeshConvertAccuracyType),
    ('Medium', adsk.fusion.MeshConvertAccuracyTypes.MediumMeshConvertAccuracyType),
    ('High', adsk.fusion.MeshConvertAccuracyTypes.HighMeshConvertAccuracyType),
    ('Precise', adsk.fusion.MeshConvertAccuracyTypes.PreciseMeshConvertAccuracyType),
]

REPAIR_MODES = [
    ('None', None),
    ('Stitch and remove artifacts', adsk.fusion.MeshRepairTypes.StitchAndRemoveMeshRepairType),
    ('Close holes', adsk.fusion.MeshRepairTypes.CloseHolesMeshRepairType),
    ('One-touch fix', adsk.fusion.MeshRepairTypes.OneTouchFixMeshRepairType),
    ('Wrap', adsk.fusion.MeshRepairTypes.WrapMeshRepairType),
    ('Rebuild', adsk.fusion.MeshRepairTypes.RebuildMeshRepairType),
]

# key, label, enum
METHODS = [
    ('form', 'Organic Form (T-Spline, built-in)', None),
    ('custom_organic', 'Organic NURBS (built-in)', None),
    ('prismatic', 'Prismatic (planes, cylinders, cones)',
     adsk.fusion.MeshConvertMethodTypes.PrismaticMeshConvertMethodType),
    ('organic', 'Organic (Form / T-Spline)',
     adsk.fusion.MeshConvertMethodTypes.OrganicMeshConvertMethodType),
    ('faceted', 'Faceted (one face per triangle)',
     adsk.fusion.MeshConvertMethodTypes.FacetedMeshConvertMethodType),
]
METHOD_BY_KEY = {k: (lab, en) for k, lab, en in METHODS}

RESOLUTIONS = ['Auto from max error', 'Low', 'Medium', 'High', 'Precise', 'By face count']


def probe_methods(comp, mesh_body):
    """Which conversion methods this licence actually allows.

    Setting meshConvertMethodType raises immediately when the required design
    extension is missing, so this costs nothing and creates no feature.
    """
    out = {}
    for key, label, enum in METHODS:
        if enum is None:                 # ours; nothing to ask the licence about
            ok = bool(python_cmd()) and os.path.exists(ORGANIC_ENGINE)
            out[key] = (ok, '' if ok else
                        'needs a Python with numpy, scipy and trimesh; see '
                        'python_path.txt')
            continue
        try:
            ci = comp.features.meshConvertFeatures.createInput([mesh_body])
            ci.meshConvertMethodType = enum
            out[key] = (True, '')
        except Exception:
            out[key] = (False, short_error(traceback.format_exc()))
    return out


def auto_tier(max_error_mm, diag_mm):
    if diag_mm <= 0:
        return 1
    r = max_error_mm / diag_mm
    if r > 1.0e-2:
        return 0
    if r > 3.0e-3:
        return 1
    if r > 1.0e-3:
        return 2
    return 3


# --------------------------------------------------------------------------- #
# mesh geometry - stdlib only
# --------------------------------------------------------------------------- #

def read_triangles(mesh_body):
    """(coords_cm, tri_idx): flat [x,y,z,...] and flat index triples."""
    pm = mesh_body.mesh
    if pm and pm.triangleCount > 0 and pm.quadCount == 0 and pm.polygonCount == 0:
        return list(pm.nodeCoordinatesAsDouble), list(pm.triangleNodeIndices)
    dm = mesh_body.displayMesh
    return list(dm.nodeCoordinatesAsDouble), list(dm.nodeIndices)


def weld(coords, idx, tol_cm=WELD_TOL_CM):
    """Merge coincident vertices so edge adjacency becomes real."""
    inv = 1.0 / tol_cm
    lut = {}
    remap = [0] * (len(coords) // 3)
    out = []
    for i in range(len(coords) // 3):
        x = coords[3 * i]
        y = coords[3 * i + 1]
        z = coords[3 * i + 2]
        key = (int(round(x * inv)), int(round(y * inv)), int(round(z * inv)))
        j = lut.get(key)
        if j is None:
            j = len(out) // 3
            lut[key] = j
            out.extend((x, y, z))
        remap[i] = j
    return out, [remap[i] for i in idx]


def boundary_loops(idx):
    """Directed half-edge walk. Returns vertex loops bounding open holes."""
    he = set()
    for t in range(0, len(idx), 3):
        a, b, c = idx[t], idx[t + 1], idx[t + 2]
        if a == b or b == c or a == c:
            continue
        he.add((a, b))
        he.add((b, c))
        he.add((c, a))

    bnd = [e for e in he if (e[1], e[0]) not in he]
    if not bnd:
        return []

    nxt = {}
    for a, b in bnd:
        nxt.setdefault(a, []).append(b)

    used = set()
    loops = []
    limit = len(bnd) + 2
    for e0 in bnd:
        if e0 in used:
            continue
        a0, b0 = e0
        used.add(e0)
        loop = [a0]
        cur = b0
        guard = 0
        while cur != a0 and guard < limit:
            loop.append(cur)
            nb = None
            for cand in nxt.get(cur, ()):
                if (cur, cand) not in used:
                    nb = cand
                    break
            if nb is None:
                break
            used.add((cur, nb))
            cur = nb
            guard += 1
        if cur == a0 and len(loop) >= 3:
            loops.append(loop)
    return loops


SPAN_SAMPLE_CAP = 600


def loop_span_mm(coords, loop):
    """Gap size in mm: the widest straight-line distance across the hole.

    Not the bounding-box diagonal - that over-reports a round hole by a factor of
    sqrt(2), which would make the gap limit mean something different for every hole
    shape. Long loops are strided down to SPAN_SAMPLE_CAP points; the answer is then
    a slight under-estimate rather than a wrong-by-construction one.
    """
    pts = loop
    if len(pts) > SPAN_SAMPLE_CAP:
        stride = len(pts) // SPAN_SAMPLE_CAP + 1
        pts = pts[::stride]
    xyz = [(coords[3 * v], coords[3 * v + 1], coords[3 * v + 2]) for v in pts]
    best = 0.0
    for i in range(len(xyz) - 1):
        ax, ay, az = xyz[i]
        for j in range(i + 1, len(xyz)):
            bx, by, bz = xyz[j]
            d = (ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2
            if d > best:
                best = d
    return math.sqrt(best) * CM


def fill_loops(coords, idx, loops, max_gap_mm, fill_all):
    """Fan-fill every loop at or under max_gap_mm. The patch winding is reversed
    against the boundary half-edge so its normal agrees with the surrounding shell."""
    out_c = list(coords)
    out_i = list(idx)
    filled, left = [], []
    for loop in loops:
        size = loop_span_mm(out_c, loop)
        if not fill_all and size > max_gap_mm:
            left.append(size)
            continue
        n = len(loop)
        cx = sum(out_c[3 * v] for v in loop) / n
        cy = sum(out_c[3 * v + 1] for v in loop) / n
        cz = sum(out_c[3 * v + 2] for v in loop) / n
        ci = len(out_c) // 3
        out_c.extend((cx, cy, cz))
        for k in range(n):
            a = loop[k]
            b = loop[(k + 1) % n]
            out_i.extend((b, a, ci))
        filled.append(size)
    return out_c, out_i, filled, left


def stats_mm(dists_cm):
    if not dists_cm:
        return None
    vals = sorted(abs(d) * CM for d in dists_cm)
    n = len(vals)
    return {
        'n': n,
        'max': vals[-1],
        'mean': sum(vals) / n,
        'rms': math.sqrt(sum(v * v for v in vals) / n),
        'p95': vals[min(n - 1, int(0.95 * n))],
        'p99': vals[min(n - 1, int(0.99 * n))],
    }


# --------------------------------------------------------------------------- #
# base-feature plumbing
# --------------------------------------------------------------------------- #

def is_parametric(design):
    return design.designType == adsk.fusion.DesignTypes.ParametricDesignType


def add_mesh_body(comp, parametric, coords, idx, name):
    """Create a mesh body from raw triangle data.

    Returns (body, owning BaseFeature, actual_name). The handle is re-fetched
    because finishEdit can invalidate it, and the name is read back because Fusion
    appends " (1)" when the name is already taken - searching for the name we asked
    for would then hand back somebody else's body.
    """
    base = None
    if parametric:
        base = comp.features.baseFeatures.add()
        base.startEdit()
    actual = name
    try:
        body = comp.meshBodies.addByTriangleMeshData(coords, idx, [], [])
        if body:
            body.name = name
            actual = body.name
    finally:
        if base:
            base.finishEdit()
    if body is None:
        return None, base, actual
    for i in range(comp.meshBodies.count - 1, -1, -1):
        b = comp.meshBodies.item(i)
        if b.name == actual:
            return b, base, actual
    return body, base, actual


def newest_brep(comp):
    return comp.bRepBodies.item(comp.bRepBodies.count - 1) if comp.bRepBodies.count else None


# --------------------------------------------------------------------------- #
# deviation measurement
# --------------------------------------------------------------------------- #

TMP_NAME = '__MeshToPart_tmp__'


def measure_deviation(comp, parametric, src_mesh_body, brep_body, max_error_mm):
    """Two-sided deviation between the converted body and the untouched original mesh.

    PolygonMesh.compareWith gives the signed closest-point distance per node natively,
    so the converted body is tessellated, injected as a throw-away mesh body, compared
    both ways, and the throw-away body removed again.
    """
    tmp_base = None
    try:
        calc = brep_body.meshManager.createMeshCalculator()
        calc.surfaceTolerance = max(max_error_mm * MM / 5.0, 1.0e-5)
        tm = calc.calculate()
        if tm is None or tm.nodeCount == 0:
            return None
        tmp, tmp_base, tmp_name = add_mesh_body(comp, parametric,
                                                list(tm.nodeCoordinatesAsDouble),
                                                list(tm.nodeIndices), TMP_NAME)
        if tmp is None:
            return None
        return {
            'out': stats_mm(tmp.mesh.compareWith(src_mesh_body.mesh, None, None)),
            'back': stats_mm(src_mesh_body.mesh.compareWith(tmp.mesh, None, None)),
        }
    except Exception:
        log('measure_deviation failed\n' + traceback.format_exc())
        return None
    finally:
        try:
            if tmp_base:
                tmp_base.deleteMe()
            for i in range(comp.meshBodies.count - 1, -1, -1):
                b = comp.meshBodies.item(i)
                if b.name == TMP_NAME or b.name.startswith(TMP_NAME):
                    b.deleteMe()
        except Exception:
            log('temp cleanup failed\n' + traceback.format_exc())


# --------------------------------------------------------------------------- #
# the escalation ladder
# --------------------------------------------------------------------------- #

def build_ladder(o, diag_mm):
    """Ordered attempts. Each rung says what to change before re-converting.

    Every method gets a ladder, because the max-error budget has to mean something
    whichever conversion the licence allows:
      organic    - raise the accuracy tier
      prismatic  - tighter face groups, then the accurate fitter
      faceted    - drop the reduction, which is the only thing that lost accuracy
    """
    n = max(1, o['max_attempts']) if o['escalate'] else 1
    key = o['method_key']
    if key == 'form':
        n = max(n, 3)          # the step-down fallbacks are not optional
    rungs = []

    if key == 'form':
        # For a Form body the binding constraint is validity, not the budget:
        # Fusion refuses a self-intersecting T-Spline outright and returns no body.
        # A finer cage follows the scan better but is likelier to self-intersect,
        # so this ladder steps DOWN until Fusion accepts the cage.
        start = int(o['cage'])
        if o.get('cage_auto') and diag_mm > 0:
            start = int(max(14, min(28, round(0.06 * diag_mm /
                                              max(o['max_error_mm'], 1e-9)))))
        seen = []
        for f in (1.0, 0.75, 0.55, 0.4):
            c = max(10, int(round(start * f)))
            if c not in seen:
                seen.append(c)
        return [{'label': 'cage %d' % c, 'cage': c} for c in seen[:max(1, n)]]

    if key == 'custom_organic':
        # accuracy is set by how finely the cage samples the mesh, so that is what
        # the ladder climbs
        start = int(o['cage'])
        if o.get('cage_auto') and diag_mm > 0:
            # calibrated on this build: p99 deviation runs at roughly
            # 0.03 * bbox_diagonal / cage for a scan-like mesh, better on clean
            # geometry. Invert that for a starting point, then let the ladder
            # correct it against the measurement.
            # measured on this build: p99 deviation runs at about
            # 0.06 * bbox_diagonal / cage on a scan-like mesh. The floor is 20
            # because below that the surface is too coarse to be worth building
            # whatever the budget says, and the ladder corrects upward from here.
            want = 0.06 * diag_mm / max(o['max_error_mm'], 1.0e-9)
            start = int(max(20, min(48, round(want))))
        for k in range(n):
            c = int(round(start * (1.35 ** k)))
            rungs.append({'label': 'cage %d' % c, 'cage': c})
        return rungs

    if key == 'organic':
        if o['resolution'] == 'By face count':
            return [{'label': 'faces=%d' % o['face_count']}]
        start = (auto_tier(o['max_error_mm'], diag_mm)
                 if o['resolution'] == 'Auto from max error'
                 else [t[0] for t in ACC_TIERS].index(o['resolution']))
        for t in range(start, min(len(ACC_TIERS), start + n)):
            rungs.append({'label': ACC_TIERS[t][0], 'tier': t})

    elif key == 'prismatic':
        angle = o['feature_angle_rad']
        btol = min(0.01, max(1.0e-5, o['max_error_mm'] / max(diag_mm, 1.0e-9)))
        cand = [
            {'label': 'groups %.0f deg' % math.degrees(angle),
             'groups': ('fast', angle, o['min_group'])},
            {'label': 'groups %.0f deg' % math.degrees(angle * 0.6),
             'groups': ('fast', angle * 0.6, o['min_group'] * 0.5)},
            {'label': 'accurate groups tol %.1e' % btol,
             'groups': ('accurate', None, btol)},
            {'label': 'no reduction', 'drop_reduce': True,
             'groups': ('accurate', None, btol)},
        ]
        rungs = cand[:n] if o['face_groups'] else [{'label': 'as imported'}]

    else:                                   # faceted
        cand = [{'label': 'as prepared'}, {'label': 'no reduction', 'drop_reduce': True}]
        rungs = cand[:n]

    return rungs or [{'label': 'default'}]


# --------------------------------------------------------------------------- #
# the pipeline
# --------------------------------------------------------------------------- #

def run_pipeline(app, src_body, o, progress=None):
    design = adsk.fusion.Design.cast(app.activeProduct)
    parametric = is_parametric(design)
    comp = src_body.parentComponent

    # A solid can be offset too, and every line below here assumes a mesh body
    # (isClosed, displayMesh, meshConvertFeatures). Tessellate once, right here,
    # so exactly one place in the pipeline knows the difference exists.
    if adsk.fusion.BRepBody.cast(src_body):
        t_coords, t_idx = read_brep_triangles(src_body)
        tess, _tb, _tn = add_mesh_body(comp, parametric, t_coords, t_idx,
                                       '%s [work tess %s]'
                                       % (src_body.name, time.strftime('%H%M%S')))
        if not tess:
            raise RuntimeError('Could not tessellate the solid body.')
        src_body = tess

    rep = {'source': src_body.name, 'warnings': [], 'attempts': [], 'version': VERSION}

    def step(msg):
        log('  ' + msg)
        if progress:
            try:
                progress.message = msg
                adsk.doEvents()
            except Exception:
                pass

    bb = src_body.boundingBox
    ext = [(bb.maxPoint.x - bb.minPoint.x) * CM,
           (bb.maxPoint.y - bb.minPoint.y) * CM,
           (bb.maxPoint.z - bb.minPoint.z) * CM]
    diag_mm = math.sqrt(sum(e * e for e in ext))
    rep['extent_mm'] = ext
    rep['diag_mm'] = diag_mm
    rep['closed_in'] = src_body.isClosed
    rep['oriented_in'] = src_body.isOriented

    # what this licence actually allows, before anything is built
    o = dict(o)
    if o['method_key'] == 'form' and not parametric:
        raise RuntimeError(
            'A Form (T-Spline) body needs a parametric design, and this one has '
            'design history captured off. Turn it on (right-click the document '
            'name in the browser, "Capture Design History") or pick Organic NURBS, '
            'which works in a direct design.')
    avail = probe_methods(comp, src_body)
    rep['methods_available'] = sorted(k for k, v in avail.items() if v[0])
    if not avail.get(o['method_key'], (False, 'unknown method'))[0]:
        why = avail.get(o['method_key'], (False, 'unknown method'))[1]
        if o['fallback_faceted'] and avail.get('faceted', (False, ''))[0]:
            rep['warnings'].append('%s is not available on this licence (%s). Using '
                                   'Faceted instead.'
                                   % (o['method_label'].split(' (')[0], why))
            o['method_key'] = 'faceted'
            o['method_label'] = METHOD_BY_KEY['faceted'][0]
            rep['fell_back'] = True
        else:
            raise RuntimeError('%s is not available on this licence.\n%s'
                               % (o['method_label'], why))

    # 1 - read, weld, find the gaps
    step('Reading mesh')
    coords, idx = read_triangles(src_body)
    rep['tri_in'] = len(idx) // 3
    step('Welding %d triangles' % rep['tri_in'])
    wcoords, widx = weld(coords, idx)
    rep['nodes_in'] = len(wcoords) // 3

    step('Finding open boundaries')
    loops = boundary_loops(widx)
    rep['loops'] = sorted((loop_span_mm(wcoords, lp) for lp in loops), reverse=True)

    work_coords, work_idx = wcoords, widx
    rep['filled'], rep['left_open'] = [], []
    if loops and o['fill_gaps']:
        step('Filling open boundaries up to %.3f mm' % o['max_gap_mm'])
        work_coords, work_idx, filled, left = fill_loops(
            wcoords, widx, loops, o['max_gap_mm'], o['fill_all'])
        rep['filled'] = sorted(filled, reverse=True)
        rep['left_open'] = sorted(left, reverse=True)

    # 2 - working copy, so the original is never touched
    step('Creating working mesh copy')
    # unique per run: a leftover "[work]" from an earlier failed run would
    # otherwise be what refetch() finds
    work_name = '%s [work %s]' % (src_body.name, time.strftime('%H%M%S'))
    work, work_base, work_name = add_mesh_body(comp, parametric, work_coords,
                                               work_idx, work_name)
    if not work:
        raise RuntimeError('Could not create the working mesh body.')
    rep['work_body'] = work_name

    def refetch():
        for i in range(comp.meshBodies.count):
            b = comp.meshBodies.item(i)
            if b.name == work_name:
                return b
        return None

    def add_feature(features, build):
        """Add a mesh feature. A failure to build the input is a real error and is
        raised straight away; only a rejected add() is retried inside the base
        feature that owns the mesh body."""
        inp = build(None)
        try:
            return features.add(inp)
        except Exception:
            if not work_base:
                raise
            log('direct add rejected, retrying inside the base feature\n'
                + traceback.format_exc())
            work_base.startEdit()
            try:
                return features.add(build(work_base))
            finally:
                work_base.finishEdit()

    # 3 - optional Fusion repair
    if o['repair_type'] is not None:
        step('Mesh repair: %s' % o['repair_label'])

        def build_repair(base):
            ri = comp.features.meshRepairFeatures.createInput(refetch() or work)
            ri.meshRepairType = o['repair_type']
            if o['repair_type'] == adsk.fusion.MeshRepairTypes.RebuildMeshRepairType:
                ri.meshRepairRebuildType = adsk.fusion.MeshRepairRebuildTypes \
                    .PreserveSharpEdgesMeshRepairRebuildType
                ri.density = adsk.core.ValueInput.createByReal(float(o['rebuild_density']))
            if base:
                ri.targetBaseFeature = base
            return ri

        try:
            add_feature(comp.features.meshRepairFeatures, build_repair)
            work = refetch() or work
        except Exception:
            rep['warnings'].append('Mesh repair (%s) failed: %s'
                                   % (o['repair_label'], short_error(traceback.format_exc())))
            log(traceback.format_exc())

    # 4 - optional reduction, the one place Fusion takes a literal mm budget
    tri_now = work.mesh.triangleCount if work.mesh else rep['tri_in']
    rep['tri_prepared'] = tri_now
    reduce_feat = None
    if o['reduce'] and tri_now > o['reduce_over']:
        dev_cm = o['reduce_share'] * o['max_error_mm'] * MM
        step('Reducing %d triangles at max deviation %.4f mm' % (tri_now, dev_cm * CM))

        def build_reduce(base):
            di = comp.features.meshReduceFeatures.createInput(refetch() or work)
            di.meshReduceMethodType = adsk.fusion.MeshReduceMethodTypes.AdaptiveReduceType
            di.meshReduceTargetType = adsk.fusion.MeshReduceTargetTypes \
                .MaximumDeviationMeshReduceTargetType
            di.maximumDeviation = adsk.core.ValueInput.createByReal(dev_cm)
            if base:
                di.targetBaseFeature = base
            return di

        try:
            reduce_feat = add_feature(comp.features.meshReduceFeatures, build_reduce)
            work = refetch() or work
            rep['tri_reduced'] = work.mesh.triangleCount if work.mesh else None
        except Exception:
            rep['warnings'].append('Mesh reduce failed: %s'
                                   % short_error(traceback.format_exc()))
            log(traceback.format_exc())

    # 4b - faceted makes one BRep face per triangle, so an unreduced scan would
    # produce a body Fusion cannot handle. Cap it rather than hang.
    def apply_faceted_cap():
        nonlocal work, reduce_feat
        w = refetch() or work
        n_tri = w.mesh.triangleCount if w.mesh else 0
        if n_tri <= o['faceted_cap']:
            return
        step('Faceted cap: reducing %d triangles to %d faces'
             % (n_tri, o['faceted_cap']))

        def build_cap(base):
            di = comp.features.meshReduceFeatures.createInput(refetch() or work)
            di.meshReduceMethodType = adsk.fusion.MeshReduceMethodTypes.AdaptiveReduceType
            di.meshReduceTargetType = adsk.fusion.MeshReduceTargetTypes \
                .FaceCountMeshReduceTargetType
            di.facecount = adsk.core.ValueInput.createByReal(float(o['faceted_cap']))
            if base:
                di.targetBaseFeature = base
            return di

        try:
            reduce_feat = add_feature(comp.features.meshReduceFeatures, build_cap)
            work = refetch() or work
            rep['faceted_capped'] = [n_tri, work.mesh.triangleCount if work.mesh else None]
            rep['warnings'].append(
                'Faceted conversion makes one face per triangle, so the mesh was '
                'reduced to %d faces first. That spends error budget, which the '
                'measurement below accounts for.' % o['faceted_cap'])
        except Exception:
            rep['warnings'].append('Faceted cap reduction failed: %s'
                                   % short_error(traceback.format_exc()))
            log(traceback.format_exc())

    if o['method_key'] == 'faceted':
        apply_faceted_cap()

    # 4b - optional offset ("wrap"). A level-set offset, so it cannot
    # self-intersect the way Offset Face and Shell do in a concavity.
    if abs(float(o.get('wrap_offset_mm', 0.0))) > 1.0e-9:
        step('Wrapping: offset %+.3f mm' % o['wrap_offset_mm'])
        try:
            w = refetch() or work
            wc, wi = read_triangles(w)
            wres = run_wrap(wc, wi,
                            {'offset_mm': float(o['wrap_offset_mm']),
                             'pitch_mm': float(o.get('wrap_pitch_mm', 0.0)),
                             'close_mm': float(o.get('wrap_close_mm', 0.0)),
                             'smooth_iters': int(o.get('wrap_smooth', 15))})
            work, work_base, work_name = add_mesh_body(
                comp, parametric, wres['coords'], wres['idx'],
                '%s [work wrap %s]' % (rep['source'], time.strftime('%H%M%S')))
            if not work:
                raise RuntimeError('Could not create the offset mesh body.')
            rep['wrap'] = wres['stats']
        except Exception:
            rep['warnings'].append('Wrap failed: %s'
                                   % short_error(traceback.format_exc()))
            log(traceback.format_exc())

    # 5 - the conversion ladder
    def apply_groups(spec):
        kind, angle, tol = spec

        def build_groups(base):
            gi = comp.features.meshGenerateFaceGroupsFeatures.createInput(refetch() or work)
            if kind == 'fast':
                gi.meshGenerateFaceGroupsMethodType = adsk.fusion \
                    .MeshGenerateFaceGroupsMethodTypes.FastGenerateFaceGroupsType
                gi.angleThreshold = adsk.core.ValueInput.createByReal(angle)
                gi.minimumFaceGroupSize = adsk.core.ValueInput.createByReal(tol)
            else:
                gi.meshGenerateFaceGroupsMethodType = adsk.fusion \
                    .MeshGenerateFaceGroupsMethodTypes.AccurateGenerateFaceGroupsType
                gi.boundaryTolerance = adsk.core.ValueInput.createByReal(tol)
            if base:
                gi.targetBaseFeature = base
            return gi

        return add_feature(comp.features.meshGenerateFaceGroupsFeatures, build_groups)

    method_enum = METHOD_BY_KEY[o['method_key']][1]
    organic = o['method_key'] == 'organic'

    def build_convert(base, rung):
        body = refetch() or work
        try:
            ci = comp.features.meshConvertFeatures.createInput([body])
        except Exception:
            oc = adsk.core.ObjectCollection.create()
            oc.add(body)
            ci = comp.features.meshConvertFeatures.createInput(oc)
        ci.meshConvertMethodType = method_enum
        ci.meshConvertOperationType = o['op_type']
        if organic:
            if 'tier' in rung:
                ci.meshConvertResolutionType = adsk.fusion.MeshConvertResolutionTypes \
                    .ByAccuracyMeshConvertResolutionType
                ci.meshConvertAccuracyType = ACC_TIERS[rung['tier']][1]
            else:
                ci.meshConvertResolutionType = adsk.fusion.MeshConvertResolutionTypes \
                    .ByFacetNumberMeshConvertResolutionType
                ci.numberOfFaces = adsk.core.ValueInput.createByReal(float(o['face_count']))
            ci.isPreprocessHoles = o['preprocess_holes']
        if base:
            ci.targetBaseFeature = base
        return ci

    ladder = build_ladder(o, diag_mm)
    rep['ladder'] = [r['label'] for r in ladder]
    forms_before = comp.features.formFeatures.count
    feat = None
    result_body = None
    dev = None
    last_error = None
    org_log = []

    def worst_mm(d):
        """Two-sided deviation. part->mesh alone is not a criterion: a body that
        covers a tenth of the scan scores zero on it. mesh->part is what catches
        geometry the conversion dropped."""
        if not d:
            return None
        a = d['out']['max'] if d.get('out') else None
        b = d['back']['max'] if d.get('back') else None
        vals = [v for v in (a, b) if v is not None]
        return max(vals) if vals else None

    def usable(body):
        """A body worth keeping even when it misses the budget."""
        if body is None:
            return False
        return body.isSolid or not src_body.isClosed

    def acceptable(body, d):
        if not usable(body):
            return False          # a closed scan must give a solid, not a patch
        if not o['verify']:
            return True
        w = worst_mm(d)
        if w is None:
            return True           # unmeasurable: do not loop forever
        return w <= o['max_error_mm']

    def try_rung(rung, label):
        """One conversion attempt. Returns (body, deviation) or (None, None) on a
        conversion failure - a refusal by the modeller is a result, not a crash,
        because the next rung exists precisely to get past it."""
        nonlocal feat, last_error
        if rung.get('groups'):
            try:
                apply_groups(rung['groups'])
                w = refetch()
                if w:
                    rep['face_groups'] = w.faceGroups.count
            except Exception:
                rep['warnings'].append('Face-group generation failed: %s'
                                       % short_error(traceback.format_exc()))
                log(traceback.format_exc())

        step('Converting: %s (%s)' % (o['method_label'].split(' (')[0], label))
        try:
            if o['method_key'] in ('custom_organic', 'form'):
                # ctrl < cage+1 makes this a least-squares FIT; ctrl == cage+1
                # makes it an interpolation of every cage point. The cage sits
                # exactly on the mesh, so interpolating it reproduces the
                # mesh's own faceting as curvature ripple -- visible as pinch
                # lines / stripes under specular shading. Sampling densely
                # (cage) and fitting loosely (ctrl) is what removes them.
                m = rung['cage'] + 1
                ctrl = max(6, min(m, int(round(o.get('ctrl_ratio', 1.0) * m))))
                params = {'cage': rung['cage'],
                          'wrap_iters': o['wrap_iters'],
                          'relax_iters': o['relax_iters'],
                          'ctrl': ctrl,
                          'fair': o['fairness'],
                          'aspect_limit': float(o.get('aspect_limit', 0.0)),
                          'aspect_iters': int(o.get('aspect_iters', 0)),
                          'aspect_gain': float(o.get('aspect_gain', 2.0))}
                w = refetch() or work
                builder = (organic_form_body if o['method_key'] == 'form'
                           else organic_body)
                body, eres, feat = builder(app, w, params, org_log)
                rep['engine'] = eres['stats']
                rep['attempts'].append({
                    'rung': label,
                    'max_mm': max(eres['stats']['surf_to_mesh_max_mm'],
                                  eres['stats']['mesh_to_surf_max_mm']),
                    'out_mm': eres['stats']['surf_to_mesh_max_mm'],
                    'back_mm': eres['stats']['mesh_to_surf_max_mm'],
                    'faces': body.faces.count, 'solid': body.isSolid,
                })
                if not o['verify']:
                    return body, None
                step('Measuring deviation against the original mesh')
                d = measure_deviation(comp, parametric, src_body, body,
                                      o['max_error_mm'])
                rep['attempts'][-1].update({
                    'max_mm': worst_mm(d),
                    'out_mm': None if not d else d['out']['max'],
                    'back_mm': None if not d or not d.get('back')
                               else d['back']['max'],
                })
                return body, d
            feat = add_feature(comp.features.meshConvertFeatures,
                               lambda base, _r=rung: build_convert(base, _r))
        except Exception:
            last_error = short_error(traceback.format_exc())
            log(traceback.format_exc())
            feat = None
            rep['attempts'].append({'rung': label, 'max_mm': None, 'faces': 0,
                                    'error': last_error})
            return None, None

        body = (feat.bodies.item(0) if (feat and feat.bodies.count)
                else newest_brep(comp))
        if body is None:
            last_error = 'the conversion reported success but produced no body'
            rep['attempts'].append({'rung': label, 'max_mm': None, 'faces': 0,
                                    'error': last_error})
            return None, None

        if not o['verify']:
            rep['attempts'].append({'rung': label, 'max_mm': None,
                                    'faces': body.faces.count,
                                    'solid': body.isSolid})
            return body, None

        step('Measuring deviation against the original mesh')
        d = measure_deviation(comp, parametric, src_body, body, o['max_error_mm'])
        rep['attempts'].append({
            'rung': label,
            'max_mm': worst_mm(d),
            'out_mm': None if not d else d['out']['max'],
            'back_mm': None if not d or not d.get('back') else d['back']['max'],
            'faces': body.faces.count,
            'solid': body.isSolid,
        })
        if d is None:
            rep['warnings'].append('Deviation could not be measured.')
        return body, d

    for n, rung in enumerate(ladder):
        if n:
            step('Retrying: %s' % rung['label'])
            try:
                if feat:
                    feat.deleteMe()
                feat = None
            except Exception:
                rep['warnings'].append('Could not roll back the previous conversion; '
                                       'keeping that result.')
                log(traceback.format_exc())
                break
            if rung.get('drop_reduce') and reduce_feat:
                try:
                    reduce_feat.deleteMe()
                    reduce_feat = None
                    rep['reduction_dropped'] = True
                except Exception:
                    log(traceback.format_exc())

        body, d = try_rung(rung, rung['label'])
        if body is None:
            result_body = None
            continue
        result_body, dev = body, d
        if o['method_key'] == 'form' and usable(body):
            # the Form ladder is ordered finest-first and steps down only because
            # Fusion may reject a cage; once a cage is accepted, a coarser one can
            # only be worse, so stop here even if the budget was missed
            break
        if acceptable(body, d):
            break

    # nothing acceptable - fall back to the one method that cannot fail on a closed
    # mesh, rather than handing back a degenerate body or nothing at all
    # Only fall back when there is nothing usable. A solid that merely misses the
    # budget is kept and reported as over budget - swapping it for a 10k-face
    # faceted body would throw away the very thing that was asked for.
    if (not usable(result_body)) and o['fallback_faceted'] \
            and o['method_key'] != 'faceted' and avail.get('faceted', (False, ''))[0]:
        rep['warnings'].append(
            '%s produced no usable body (%s). Fell back to Faceted, which is exact on '
            'a closed mesh but makes one face per triangle.'
            % (o['method_label'].split(' (')[0],
               last_error or 'the result was not a solid'))
        method_enum = METHOD_BY_KEY['faceted'][1]
        organic = False
        o['method_key'] = 'faceted'
        o['method_label'] = METHOD_BY_KEY['faceted'][0]
        rep['fell_back'] = True
        try:
            if feat:
                feat.deleteMe()
                feat = None
        except Exception:
            log(traceback.format_exc())
        apply_faceted_cap()
        for lbl in ('faceted fallback', 'faceted fallback, no reduction'):
            if lbl.endswith('no reduction'):
                try:
                    if feat:
                        feat.deleteMe()
                        feat = None
                except Exception:
                    log(traceback.format_exc())
                if reduce_feat:
                    try:
                        reduce_feat.deleteMe()
                        reduce_feat = None
                        rep['reduction_dropped'] = True
                    except Exception:
                        log(traceback.format_exc())
            body, d = try_rung({'label': lbl}, lbl)
            if body is not None:
                result_body, dev = body, d
            if acceptable(body, d):
                break
        rep['ladder'].append('faceted fallback')

    if result_body is None:
        try:
            if work_base:
                work_base.deleteMe()          # do not leave a stray prepared mesh
        except Exception:
            log(traceback.format_exc())
        raise RuntimeError(
            'No conversion succeeded. Last refusal: %s\nMethods available on this '
            'licence: %s' % (last_error,
                             ', '.join(k for k, v in avail.items() if v[0]) or 'none'))

    rep['rung_used'] = rep['attempts'][-1]['rung'] if rep['attempts'] else ladder[0]['label']
    rep['deviation'] = dev
    rep['is_form'] = comp.features.formFeatures.count > forms_before
    rep['faces_out'] = result_body.faces.count
    rep['solid_out'] = result_body.isSolid
    rep['volume_out_mm3'] = result_body.volume * 1000.0 if result_body.isSolid else None
    rep['volume_in_mm3'] = src_body.volume * 1000.0 if src_body.isClosed else None

    rep['method_label'] = o['method_label']
    rep['method_key'] = o['method_key']
    result_body.name = '%s [part]' % src_body.name
    rep['result'] = result_body.name

    if org_log:
        rep['engine_log'] = org_log[-30:]
    leftover = refetch()
    if leftover and o['method_key'] == 'custom_organic':
        try:
            leftover.isLightBulbOn = False
        except Exception:
            pass
        leftover = None
    if leftover:
        try:
            leftover.isLightBulbOn = False
        except Exception:
            pass
        rep['warnings'].append('The prepared mesh "%s" survived the conversion; it is '
                               'hidden, delete it once you are happy.' % work_name)

    if o['hide_source']:
        try:
            src_body.isLightBulbOn = False
        except Exception:
            pass

    return rep


# --------------------------------------------------------------------------- #
# report
# --------------------------------------------------------------------------- #

def format_report(rep, o, secs):
    L = []
    L.append('%s  ->  %s' % (rep['source'], rep.get('result', '(none)')))
    L.append('')
    L.append('INPUT')
    L.append('  triangles         %d' % rep['tri_in'])
    L.append('  nodes (welded)    %d' % rep['nodes_in'])
    if rep.get('extent_mm'):
        L.append('  bounding box      %.2f x %.2f x %.2f mm' % tuple(rep['extent_mm']))
    L.append('  bbox diagonal     %.3f mm' % rep['diag_mm'])
    L.append('  closed / oriented %s / %s' % (rep['closed_in'], rep['oriented_in']))
    if rep['loops']:
        top = ', '.join('%.3f' % v for v in rep['loops'][:8])
        L.append('  open boundaries   %d  (mm: %s%s)'
                 % (len(rep['loops']), top, ' ...' if len(rep['loops']) > 8 else ''))
    else:
        L.append('  open boundaries   none')

    L.append('')
    L.append('PREPARATION')
    if rep['filled']:
        L.append('  gaps filled       %d  (largest %.3f mm, limit %.3f mm)'
                 % (len(rep['filled']), max(rep['filled']), o['max_gap_mm']))
    elif rep['loops']:
        L.append('  gaps filled       0')
    if rep['left_open']:
        L.append('  left open         %d  (largest %.3f mm, above the gap limit)'
                 % (len(rep['left_open']), max(rep['left_open'])))
    if rep.get('tri_reduced') is not None:
        L.append('  reduced           %d -> %d triangles at <= %.4f mm'
                 % (rep['tri_prepared'], rep['tri_reduced'],
                    o['reduce_share'] * o['max_error_mm']))
    else:
        L.append('  reduction         not applied')
    if rep.get('faceted_capped'):
        L.append('  faceted cap       %s -> %s triangles'
                 % tuple(rep['faceted_capped']))
    if rep.get('reduction_dropped'):
        L.append('  reduction         dropped on retry to recover accuracy')
    if rep.get('face_groups') is not None and rep.get('method_key') != 'faceted':
        L.append('  face groups       %d' % rep['face_groups'])

    L.append('')
    if rep.get('wrap'):
        w = rep['wrap']
        a = w['achieved']
        L.append('')
        L.append('WRAP (level-set offset)')
        L.append('  requested         %+.3f mm' % w['requested_offset_mm'])
        L.append('  achieved          mean %+.4f   min %+.4f   max %+.4f mm'
                 % (a['mean_mm'], a['min_mm'], a['max_mm']))
        L.append('  spread            rms %.4f mm about that mean, %d points'
                 % (a['rms_err_mm'], a['points']))
        L.append('  voxel pitch       %.4f mm  (%.1f M voxels)'
                 % (w['pitch_mm'], w['voxels'] / 1.0e6))
        if w['close_mm']:
            L.append('  gap closing       %.3f mm' % w['close_mm'])
        L.append('  offset mesh       %d triangles, watertight %s'
                 % (w['out_tris'], w['watertight']))
        L.append('  NOTE the deviation below is measured against the OFFSET mesh,')
        L.append('       not the body you selected.')


    L.append('CONVERSION')
    L.append('  method            %s%s'
             % (rep.get('method_label', o['method_label']),
                '   (fell back from %s)' % o['method_label'].split(' (')[0]
                if rep.get('fell_back') else ''))
    if rep.get('methods_available') is not None:
        L.append('  licence allows    %s' % (', '.join(rep['methods_available'])
                                             or 'nothing'))
    L.append('  ladder            %s' % ' -> '.join(rep['ladder']))
    L.append('  used              %s' % rep['rung_used'])
    L.append('  faces produced    %d' % rep['faces_out'])
    L.append('  solid             %s' % rep['solid_out'])
    L.append('  editable in Form  %s' % ('yes' if rep['is_form'] else 'no'))
    if rep.get('engine'):
        e = rep['engine']
        form = e.get('output') == 'tsm'
        L.append('')
        L.append('  %s ENGINE' % ('ORGANIC FORM (T-SPLINE)' if form
                                  else 'ORGANIC NURBS'))
        L.append('    cage            %d per cube face, %d points'
                 % (e['cage'], e['cage_verts']))
        if form:
            L.append('    control mesh    %d quads, %d folded'
                     % (e.get('cage_quads', 0), e.get('cage_folded_quads', 0)))
            L.append('    limit fit       surface sits %.4f mm from the mesh at the'
                     ' cage points' % e.get('limit_at_cage_max_mm', 0.0))
            L.append('    regularisation  %.1e' % e.get('limit_regularisation', 0.0))
        else:
            L.append('    control net     %d x %d per patch, fairness %.3g'
                     % (e['control_per_side'], e['control_per_side'], e['fairness']))
            L.append('    cage lands on the mesh to %.2e mm'
                     % e.get('cage_on_mesh_max_mm', 0.0))
            if e.get('cage_folded_quads') is not None:
                L.append('    folded quads    %d of the cage (intrinsic test)'
                         % e['cage_folded_quads'])
        q = e.get('cage_quality')
        if q and q.get('after'):
            aft = q['after']
            bef = q.get('before')
            L.append('    quad shape      aspect med %.2f  p95 %.2f  max %.2f'
                     '   skew p95 %.2f'
                     % (aft['aspect_median'], aft['aspect_p95'],
                        aft['aspect_max'], aft['skew_p95']))
            if bef:
                L.append('                    was       med %.2f  p95 %.2f  max %.2f'
                         '   skew p95 %.2f  (limit %.2f)'
                         % (bef['aspect_median'], bef['aspect_p95'],
                            bef['aspect_max'], bef['skew_p95'], aft['limit']))
            L.append('    %d of %d quads still exceed the aspect limit'
                     % (aft['over_limit'], aft['quads']))
        if e.get('signed_volume_mm3'):
            L.append('    closed volume   %.1f vs mesh %.1f mm3  (%+.3f %%)'
                     % (e['signed_volume_mm3'], e['mesh_volume_mm3'],
                        100.0 * (e['signed_volume_mm3'] - e['mesh_volume_mm3'])
                        / e['mesh_volume_mm3']))
        L.append('    engine check    surf->mesh max %.4f p99 %.4f rms %.4f mm'
                 % (e['surf_to_mesh_max_mm'], e['surf_to_mesh_p99_mm'],
                    e['surf_to_mesh_rms_mm']))
        L.append('    within 1 mm     %.2f %% of the surface' % e['within_1_pct'])
        L.append('    engine time     %.1f s' % e['secs'])

    if len(rep['attempts']) > 1:
        L.append('')
        L.append('  attempt                        part->mesh   mesh->part   faces   solid')
        for a in rep['attempts']:
            rung = a['rung'] if len(a['rung']) <= 28 else a['rung'][:27] + '~'
            if a.get('error'):
                L.append('    %-28s refused: %s' % (rung, a['error'][:44]))
                continue
            f = lambda v: '-' if v is None else '%.4f' % v
            L.append('    %-28s %-12s %-12s %-7d %s'
                     % (rung, f(a.get('out_mm')), f(a.get('back_mm')),
                        a['faces'], a.get('solid')))

    L.append('')
    L.append('MEASURED DEVIATION vs the original mesh')
    dev = rep.get('deviation')
    if not dev:
        L.append('  not measured')
    else:
        for key, title in (('out', 'part -> mesh'), ('back', 'mesh -> part')):
            s = dev.get(key)
            if not s:
                continue
            L.append('  %-12s  max %.4f   p99 %.4f   p95 %.4f   rms %.4f mm   (%d pts)'
                     % (title, s['max'], s['p99'], s['p95'], s['rms'], s['n']))
        worst = max(v for v in (dev['out']['max'] if dev.get('out') else None,
                                dev['back']['max'] if dev.get('back') else None)
                    if v is not None)
        L.append('')
        L.append('  worst of both directions: %.4f mm' % worst)
        if rep['closed_in'] and not rep['solid_out']:
            L.append('  WARNING: the source mesh was closed but the result is not a')
            L.append('  solid, so the conversion dropped geometry.')
        if worst <= o['max_error_mm']:
            L.append('  RESULT: within the %.4f mm budget.' % o['max_error_mm'])
        else:
            L.append('  RESULT: OVER the %.4f mm budget by %.4f mm.'
                     % (o['max_error_mm'], worst - o['max_error_mm']))
            L.append('  Raise Max attempts, lower the reduction share, allow the')
            L.append('  Faceted fallback, or accept a larger budget.')

    if rep['warnings']:
        L.append('')
        L.append('NOTES')
        for w in rep['warnings']:
            L.append('  - ' + w)

    L.append('')
    L.append('elapsed %.1f s   MeshToPart %s' % (secs, rep.get('version', VERSION)))
    return '\n'.join(L)


# --------------------------------------------------------------------------- #
# settings
# --------------------------------------------------------------------------- #

DEFAULTS = {
    'maxError': 0.1, 'escalate': True, 'maxAttempts': 3,
    'fillGaps': True, 'maxGap': 1.0, 'fillAll': False,
    'repair': 0, 'rebuildDensity': 128,
    'reduce': True, 'reduceOver': 150000, 'reduceShare': 0.4, 'facetedCap': 20000,
    'methodKey': 'prismatic', 'resolution': 0, 'faceCount': 5000,
    'preprocessHoles': True, 'faceGroups': True, 'featureAngle': 25.0, 'minGroup': 0.02,
    'opType': 0, 'verify': True, 'hideSource': True, 'showReport': True,
    'cage': 24, 'fairness': 0.01, 'wrapIters': 50, 'relaxIters': 70,
    'ctrlRatio': 0.5, 'wrapOffset': 0.0, 'wrapClose': 0.0,
    'aspectLimit': 3.0, 'aspectPasses': 80,
    'cageAuto': True,
    'fallbackFaceted': True,
}


def load_settings():
    s = dict(DEFAULTS)
    try:
        with open(SETTINGS_PATH, 'r', encoding='utf-8') as fh:
            s.update(json.load(fh))
    except Exception:
        pass
    return s


def save_settings(d):
    try:
        with open(SETTINGS_PATH, 'w', encoding='utf-8') as fh:
            json.dump(d, fh, indent=2)
    except Exception:
        pass


# --------------------------------------------------------------------------- #
# command dialog
# --------------------------------------------------------------------------- #

def build_inputs(inputs, s):
    sel = inputs.addSelectionInput('selMesh', 'Body',
                                   'Select a mesh body, or a solid to offset')
    sel.addSelectionFilter('MeshBodies')
    sel.addSelectionFilter('SolidBodies')
    sel.setSelectionLimits(1, 0)

    info = inputs.addTextBoxCommandInput('info', 'Mesh', '', 7, True)
    info.isFullWidth = True
    inputs.addBoolValueInput('btnAnalyze', 'Analyse gaps', False, '', False)

    g = inputs.addGroupCommandInput('grpAcc', 'Accuracy').children
    g.addValueInput('maxError', 'Maximum error', 'mm',
                    adsk.core.ValueInput.createByReal(s['maxError'] * MM))
    g.addBoolValueInput('escalate', 'Keep retrying until the budget is met', True, '',
                        s['escalate'])
    g.addIntegerSpinnerCommandInput('maxAttempts', 'Max attempts', 1, 4, 1,
                                    int(s['maxAttempts']))

    gp = inputs.addGroupCommandInput('grpPrep', 'Mesh preparation')
    gp.isExpanded = True
    p = gp.children
    p.addBoolValueInput('fillGaps', 'Fill small gaps / holes', True, '', s['fillGaps'])
    p.addValueInput('maxGap', 'Ignore gaps larger than', 'mm',
                    adsk.core.ValueInput.createByReal(s['maxGap'] * MM))
    p.addBoolValueInput('fillAll', 'Fill every hole regardless of size', True, '',
                        s['fillAll'])
    dd = p.addDropDownCommandInput('repair', 'Fusion repair',
                                   adsk.core.DropDownStyles.TextListDropDownStyle)
    for i, (label, _) in enumerate(REPAIR_MODES):
        dd.listItems.add(label, i == int(s['repair']))
    p.addIntegerSpinnerCommandInput('rebuildDensity', 'Rebuild density', 8, 256, 8,
                                    int(s['rebuildDensity']))
    p.addBoolValueInput('reduce', 'Reduce dense meshes', True, '', s['reduce'])
    p.addIntegerSpinnerCommandInput('reduceOver', 'Only above N triangles', 1000, 5000000,
                                    10000, int(s['reduceOver']))
    p.addFloatSpinnerCommandInput('reduceShare', 'Reduction share of the budget', '',
                                  0.05, 0.9, 0.05, float(s['reduceShare']))

    gc = inputs.addGroupCommandInput('grpConv', 'Conversion')
    gc.isExpanded = True
    c = gc.children
    # populated for real once a mesh body is selected and the licence is probed
    c.addDropDownCommandInput('method', 'Method',
                              adsk.core.DropDownStyles.TextListDropDownStyle)
    dr = c.addDropDownCommandInput('resolution', 'Organic resolution',
                                   adsk.core.DropDownStyles.TextListDropDownStyle)
    for i, label in enumerate(RESOLUTIONS):
        dr.listItems.add(label, i == int(s['resolution']))
    c.addIntegerSpinnerCommandInput('faceCount', 'Face count', 50, 200000, 100,
                                    int(s['faceCount']))
    c.addIntegerSpinnerCommandInput('facetedCap', 'Faceted: cap faces at', 500, 500000,
                                    1000, int(s['facetedCap']))
    c.addBoolValueInput('cageAuto', 'Organic: cage from the error budget',
                        True, '', s['cageAuto'])
    c.addIntegerSpinnerCommandInput('cage', 'Organic: cage per face', 8, 72, 4,
                                    int(s['cage']))
    c.addFloatSpinnerCommandInput(
        'ctrlRatio', 'Organic: control net / cage (NURBS only)', '',
        0.2, 1.0, 0.05, float(s['ctrlRatio']))
    c.addFloatSpinnerCommandInput(
        'aspectLimit', 'Organic: max quad aspect (0 = off)', '',
        0.0, 20.0, 0.5, float(s['aspectLimit']))
    c.addIntegerSpinnerCommandInput('aspectPasses', 'Organic: quad shape passes',
                                    0, 400, 10, int(s['aspectPasses']))
    c.addValueInput('wrapOffset', 'Wrap: offset distance (0 = off)', 'mm',
                    adsk.core.ValueInput.createByReal(float(s['wrapOffset']) * MM))
    c.addValueInput('wrapClose', 'Wrap: close gaps up to', 'mm',
                    adsk.core.ValueInput.createByReal(float(s['wrapClose']) * MM))
    c.addFloatSpinnerCommandInput('fairness', 'Organic: fairness (NURBS only)', '',
                                  0.0, 1.0, 0.005, float(s['fairness']))
    c.addIntegerSpinnerCommandInput('relaxIters', 'Organic: relax passes',
                                    10, 300, 10, int(s['relaxIters']))
    c.addBoolValueInput('preprocessHoles', 'Smooth remaining open hole boundaries', True,
                        '', s['preprocessHoles'])
    c.addBoolValueInput('faceGroups', 'Generate face groups (prismatic)', True, '',
                        s['faceGroups'])
    c.addValueInput('featureAngle', 'Feature angle', 'deg',
                    adsk.core.ValueInput.createByReal(math.radians(s['featureAngle'])))
    c.addFloatSpinnerCommandInput('minGroup', 'Smallest face group (area fraction)', '',
                                  0.0, 1.0, 0.005, float(s['minGroup']))
    do = c.addDropDownCommandInput('opType', 'Output',
                                   adsk.core.DropDownStyles.TextListDropDownStyle)
    do.listItems.add('Parametric feature', int(s['opType']) == 0)
    do.listItems.add('Base feature (direct)', int(s['opType']) == 1)

    ga = inputs.addGroupCommandInput('grpPost', 'After conversion')
    ga.isExpanded = True
    a = ga.children
    a.addBoolValueInput('verify', 'Measure the real deviation', True, '', s['verify'])
    a.addBoolValueInput('fallbackFaceted', 'Fall back to Faceted if the method fails',
                        True, '', s['fallbackFaceted'])
    a.addBoolValueInput('hideSource', 'Hide the original mesh', True, '', s['hideSource'])
    a.addBoolValueInput('showReport', 'Show the report when finished', True, '',
                        s['showReport'])


# NO CROSS-EVENT MODULE STATE. This module is re-imported fresh for every command
# event, so anything stashed in a global while the dialog is being built is gone by
# the time the user presses OK. Everything the handlers need is read back out of the
# dialog itself. _STATE is the one exception: a module object parked in sys.modules,
# which reload does not touch, and it is only used by the dev bridge.
_STATE = sys.modules.setdefault('_mtp_state', types.ModuleType('_mtp_state'))


def populate_methods(inputs, comp, mesh_body, preferred):
    """Probe the licence and rebuild the method dropdown from what actually works."""
    avail = probe_methods(comp, mesh_body)
    dd = inputs.itemById('method')
    dd.listItems.clear()
    keys = [k for k, _l, _e in METHODS if avail[k][0]]
    if not keys:                          # should not happen; keep the dialog usable
        keys = [k for k, _l, _e in METHODS]
    for i, key in enumerate(keys):
        want = (key == preferred) if preferred in keys else (i == 0)
        dd.listItems.add(METHOD_BY_KEY[key][0], want)
    if dd.selectedItem is None and dd.listItems.count:
        dd.listItems.item(0).isSelected = True
    return avail


def probe_note(bodies):
    """Kept for the dev bridge; the dialog no longer shows it. A method the
    licence does not allow is simply absent from the list rather than explained
    with a RuntimeError the user can do nothing about."""
    if not bodies:
        return 'select a mesh body to check which conversions this licence allows'
    try:
        avail = probe_methods(bodies[0].parentComponent, bodies[0])
    except Exception:
        return ''
    blocked = ['not available: %s - %s' % (METHOD_BY_KEY[k][0], v[1])
               for k, v in avail.items() if not v[0]]
    return '<br>'.join(blocked)


def selected_method_key(inputs):
    """Resolved from the item LABEL, never from a remembered index-to-key list."""
    dd = inputs.itemById('method')
    label = dd.selectedItem.name if dd.selectedItem else ''
    for key, lab, _enum in METHODS:
        if lab == label:
            return key
    return 'faceted'


def refresh_visibility(inputs):
    key = selected_method_key(inputs)
    organic = key == 'organic'
    prismatic = key == 'prismatic'
    res = inputs.itemById('resolution')
    res.isVisible = organic
    by_count = organic and res.selectedItem and res.selectedItem.index == 5
    inputs.itemById('faceCount').isVisible = by_count
    inputs.itemById('facetedCap').isVisible = key == 'faceted'
    custom = key in ('custom_organic', 'form')
    for cid in ('cageAuto', 'fairness', 'relaxIters', 'ctrlRatio',
                'aspectLimit', 'aspectPasses'):
        inputs.itemById(cid).isVisible = custom
    inputs.itemById('cage').isVisible = custom and not inputs.itemById('cageAuto').value
    inputs.itemById('preprocessHoles').isVisible = organic
    inputs.itemById('escalate').isVisible = not by_count
    inputs.itemById('maxAttempts').isVisible = (not by_count
                                                and inputs.itemById('escalate').value)
    inputs.itemById('faceGroups').isVisible = prismatic
    fg = prismatic and inputs.itemById('faceGroups').value
    inputs.itemById('featureAngle').isVisible = fg
    inputs.itemById('minGroup').isVisible = fg

    rep_idx = inputs.itemById('repair').selectedItem.index
    inputs.itemById('rebuildDensity').isVisible = REPAIR_MODES[rep_idx][0] == 'Rebuild'
    fill = inputs.itemById('fillGaps').value
    inputs.itemById('fillAll').isVisible = fill
    inputs.itemById('maxGap').isVisible = fill and not inputs.itemById('fillAll').value
    red = inputs.itemById('reduce').value
    inputs.itemById('reduceOver').isVisible = red
    inputs.itemById('reduceShare').isVisible = red


def resolve_options(inputs):
    key = selected_method_key(inputs)
    rep_idx = inputs.itemById('repair').selectedItem.index
    op = inputs.itemById('opType').selectedItem.index
    res = inputs.itemById('resolution')
    return {
        'max_error_mm': inputs.itemById('maxError').value * CM,
        'escalate': inputs.itemById('escalate').value,
        'max_attempts': inputs.itemById('maxAttempts').value,
        'fill_gaps': inputs.itemById('fillGaps').value,
        'max_gap_mm': inputs.itemById('maxGap').value * CM,
        'fill_all': inputs.itemById('fillAll').value,
        'repair_type': REPAIR_MODES[rep_idx][1],
        'repair_label': REPAIR_MODES[rep_idx][0],
        'rebuild_density': inputs.itemById('rebuildDensity').value,
        'reduce': inputs.itemById('reduce').value,
        'reduce_over': inputs.itemById('reduceOver').value,
        'reduce_share': inputs.itemById('reduceShare').value,
        'faceted_cap': inputs.itemById('facetedCap').value,
        'cage': inputs.itemById('cage').value,
        'cage_auto': inputs.itemById('cageAuto').value,
        'ctrl_ratio': inputs.itemById('ctrlRatio').value,
        'aspect_limit': inputs.itemById('aspectLimit').value,
        'aspect_iters': inputs.itemById('aspectPasses').value,
        'wrap_offset_mm': inputs.itemById('wrapOffset').value * CM,
        'wrap_close_mm': inputs.itemById('wrapClose').value * CM,
        'fairness': inputs.itemById('fairness').value,
        'wrap_iters': 50,
        'relax_iters': inputs.itemById('relaxIters').value,
        'fallback_faceted': inputs.itemById('fallbackFaceted').value,
        'method_key': key,
        'method_label': METHOD_BY_KEY[key][0],
        'resolution': RESOLUTIONS[res.selectedItem.index if res.selectedItem else 0],
        'face_count': inputs.itemById('faceCount').value,
        'preprocess_holes': inputs.itemById('preprocessHoles').value,
        'face_groups': inputs.itemById('faceGroups').value,
        'feature_angle_rad': inputs.itemById('featureAngle').value,
        'min_group': inputs.itemById('minGroup').value,
        'op_type': (adsk.fusion.MeshConvertOperationTypes
                    .ParametricFeatureMeshConvertOperationType if op == 0 else
                    adsk.fusion.MeshConvertOperationTypes
                    .BaseFeatureMeshConvertOperationType),
        'verify': inputs.itemById('verify').value,
        'hide_source': inputs.itemById('hideSource').value,
        'show_report': inputs.itemById('showReport').value,
    }


def settings_from(inputs):
    return {
        'maxError': inputs.itemById('maxError').value * CM,
        'escalate': inputs.itemById('escalate').value,
        'maxAttempts': inputs.itemById('maxAttempts').value,
        'fillGaps': inputs.itemById('fillGaps').value,
        'maxGap': inputs.itemById('maxGap').value * CM,
        'fillAll': inputs.itemById('fillAll').value,
        'repair': inputs.itemById('repair').selectedItem.index,
        'rebuildDensity': inputs.itemById('rebuildDensity').value,
        'reduce': inputs.itemById('reduce').value,
        'reduceOver': inputs.itemById('reduceOver').value,
        'reduceShare': inputs.itemById('reduceShare').value,
        'facetedCap': inputs.itemById('facetedCap').value,
        'cage': inputs.itemById('cage').value,
        'cageAuto': inputs.itemById('cageAuto').value,
        'ctrlRatio': inputs.itemById('ctrlRatio').value,
        'aspectLimit': inputs.itemById('aspectLimit').value,
        'aspectPasses': inputs.itemById('aspectPasses').value,
        'wrapOffset': inputs.itemById('wrapOffset').value * CM,
        'wrapClose': inputs.itemById('wrapClose').value * CM,
        'fairness': inputs.itemById('fairness').value,
        'relaxIters': inputs.itemById('relaxIters').value,
        'fallbackFaceted': inputs.itemById('fallbackFaceted').value,
        'methodKey': selected_method_key(inputs),
        'resolution': (inputs.itemById('resolution').selectedItem.index
                       if inputs.itemById('resolution').selectedItem else 0),
        'faceCount': inputs.itemById('faceCount').value,
        'preprocessHoles': inputs.itemById('preprocessHoles').value,
        'faceGroups': inputs.itemById('faceGroups').value,
        'featureAngle': math.degrees(inputs.itemById('featureAngle').value),
        'minGroup': inputs.itemById('minGroup').value,
        'opType': inputs.itemById('opType').selectedItem.index,
        'verify': inputs.itemById('verify').value,
        'hideSource': inputs.itemById('hideSource').value,
        'showReport': inputs.itemById('showReport').value,
    }


def selected_bodies(sel_input):
    out = []
    for i in range(sel_input.selectionCount):
        b = adsk.fusion.MeshBody.cast(sel_input.selection(i).entity)
        if b:
            out.append(b)
    return out


INFO_CSS = ('<style>table.m{border-collapse:collapse;font-size:11px}'
            'table.m td{padding:1px 10px 1px 0;vertical-align:top}'
            'table.m td.k{color:#666}</style>')


def _rows(pairs):
    body = ''.join('<tr><td class="k">%s</td><td>%s</td></tr>' % (k, v)
                   for k, v in pairs)
    return INFO_CSS + '<table class="m">' + body + '</table>'


def quick_info(sel_input):
    bodies = selected_bodies(sel_input)
    if not bodies:
        return _rows([('', 'Select a mesh body to convert.')])
    pairs = []
    for n, b in enumerate(bodies):
        bb = b.boundingBox
        pm = b.mesh
        if len(bodies) > 1:
            pairs.append(('body', b.name))
        else:
            pairs.append(('mesh body', b.name))
        pairs.append(('triangles', '{:,}'.format(pm.triangleCount).replace(',', ' ')))
        pairs.append(('nodes', '{:,}'.format(pm.nodeCount).replace(',', ' ')))
        pairs.append(('size', '%.2f &times; %.2f &times; %.2f mm'
                      % ((bb.maxPoint.x - bb.minPoint.x) * CM,
                         (bb.maxPoint.y - bb.minPoint.y) * CM,
                         (bb.maxPoint.z - bb.minPoint.z) * CM)))
        pairs.append(('watertight', 'yes' if b.isClosed else 'no'))
        pairs.append(('oriented', 'yes' if b.isOriented else 'no'))
        if b.isClosed:
            pairs.append(('volume', '%.1f mm&sup3;' % (b.volume * 1000.0)))
    return _rows(pairs)


def deep_info(sel_input):
    bodies = selected_bodies(sel_input)
    if not bodies:
        return _rows([('', 'Select a mesh body to analyse.')])
    pairs = []
    for b in bodies:
        coords, idx = read_triangles(b)
        wc, wi = weld(coords, idx)
        sizes = sorted((loop_span_mm(wc, lp) for lp in boundary_loops(wi)),
                       reverse=True)
        pairs.append(('mesh body', b.name))
        pairs.append(('triangles', '{:,}'.format(len(idx) // 3).replace(',', ' ')))
        pairs.append(('welded nodes', '{:,}'.format(len(wc) // 3).replace(',', ' ')))
        pairs.append(('watertight', 'yes' if b.isClosed else 'no'))
        if sizes:
            pairs.append(('open boundaries', '%d' % len(sizes)))
            pairs.append(('widest gap', '%.3f mm' % sizes[0]))
            shown = ', '.join('%.3f' % v for v in sizes[:8])
            pairs.append(('gap spans',
                          shown + (' &hellip;' if len(sizes) > 8 else '')))
            pairs.append(('', 'set the gap limit above %.3f mm to close them all'
                          % sizes[0]))
        else:
            pairs.append(('open boundaries', 'none'))
    return _rows(pairs)


def collect_mesh_bodies(comp, out):
    """Every mesh body in the assembly, at ANY depth.

    The previous version walked only `comp.occurrences`, which is the direct
    children. A mesh sitting two components down -- what an imported STL
    assembly usually looks like -- was invisible, and the add-in reported
    "no mesh body in the active document" with the mesh plainly on screen.
    `allOccurrences` is flat over the whole tree and each item already
    carries its own assembly context.
    """
    for i in range(comp.meshBodies.count):
        out.append(comp.meshBodies.item(i))
    try:
        occs = comp.allOccurrences
    except Exception:
        return
    for i in range(occs.count):
        occ = occs.item(i)
        try:
            bodies = occ.component.meshBodies
            for j in range(bodies.count):
                out.append(bodies.item(j).createForAssemblyContext(occ))
        except Exception:
            pass


def on_created(app, cmd):
    """Everything the CommandCreated handler needs, so the shell stays thin."""
    cmd.isExecutedWhenPreEmpted = False
    cmd.okButtonText = 'Convert'
    cmd.setDialogInitialSize(470, 760)

    s = load_settings()
    build_inputs(cmd.commandInputs, s)

    design = adsk.fusion.Design.cast(app.activeProduct)
    sel = cmd.commandInputs.itemById('selMesh')
    found = []
    if design:
        collect_mesh_bodies(design.rootComponent, found)
    if len(found) == 1:
        sel.addSelection(found[0])
    if found:
        populate_methods(cmd.commandInputs, found[0].parentComponent, found[0],
                         s.get('methodKey'))
    else:
        populate_methods_placeholder(cmd.commandInputs, s.get('methodKey'))
    cmd.commandInputs.itemById('info').formattedText = quick_info(sel)
    refresh_visibility(cmd.commandInputs)
    dump_dialog(cmd.commandInputs)      # so the dialog can be tested, not just eyeballed
    _STATE.live_cmd = cmd            # dev bridge only; see the note on _STATE


DIALOG_DUMP = os.path.join(ADDIN_DIR, 'devjobs', 'dialog_dump.json')


def dump_dialog(inputs):
    """Snapshot every input's id, visibility and value, plus the options they
    resolve to. Written on every dialog open; only the dev bridge reads it."""
    try:
        rows = []
        for i in range(inputs.count):
            ci = inputs.item(i)
            row = {'id': ci.id, 'type': type(ci).__name__, 'visible': ci.isVisible}
            for attr in ('value', 'selectionCount'):
                try:
                    row[attr] = getattr(ci, attr)
                except Exception:
                    pass
            try:
                if ci.selectedItem:
                    row['selected'] = ci.selectedItem.name
                    row['items'] = [ci.listItems.item(k).name
                                    for k in range(ci.listItems.count)]
            except Exception:
                pass
            rows.append(row)
        out = {'ts': time.time(), 'inputs': rows}
        try:
            out['resolved'] = {k: str(v) for k, v in resolve_options(inputs).items()}
        except Exception:
            out['resolved_error'] = short_error(traceback.format_exc())
        os.makedirs(os.path.dirname(DIALOG_DUMP), exist_ok=True)
        with open(DIALOG_DUMP, 'w', encoding='utf-8') as fh:
            json.dump(out, fh, indent=2, default=str)
    except Exception:
        log('dump_dialog failed\n' + traceback.format_exc())


def populate_methods_placeholder(inputs, preferred):
    keys = [k for k, _l, _e in METHODS]
    dd = inputs.itemById('method')
    dd.listItems.clear()
    for i, key in enumerate(keys):
        dd.listItems.add(METHOD_BY_KEY[key][0],
                         key == preferred if preferred in keys else i == 0)


def on_input_changed(app, inputs, changed_id):
    if changed_id == 'btnAnalyze':
        inputs.itemById('info').formattedText = deep_info(inputs.itemById('selMesh'))
        return
    if changed_id == 'selMesh':
        bodies = selected_bodies(inputs.itemById('selMesh'))
        if bodies:
            populate_methods(inputs, bodies[0].parentComponent, bodies[0],
                             selected_method_key(inputs))
        inputs.itemById('info').formattedText = quick_info(inputs.itemById('selMesh'))
    refresh_visibility(inputs)


def on_execute(app, ui, inputs):
    o = resolve_options(inputs)
    save_settings(settings_from(inputs))
    bodies = selected_bodies(inputs.itemById('selMesh'))
    if not bodies:
        ui.messageBox('No mesh body selected.', 'Mesh to Part')
        return

    progress = ui.createProgressDialog()
    progress.show('Mesh to Part', 'Starting', 0, len(bodies), 0)
    reports, failed = [], []
    t0 = time.time()
    for n, b in enumerate(bodies):
        progress.progressValue = n
        log('=== %s (%s) ===' % (b.name, o['method_key']))
        try:
            rep = run_pipeline(app, b, o, progress)
            reports.append(format_report(rep, o, time.time() - t0))
        except Exception:
            tb = traceback.format_exc()
            log(tb)
            failed.append(b.name)
            reports.append('%s: FAILED\n  %s\n\nFull traceback in\n  %s'
                           % (b.name, short_error(tb), LOG_PATH))
    progress.hide()

    text = ('\n\n' + '-' * 66 + '\n\n').join(reports)
    path = write_report(text)
    if o['show_report']:
        ui.messageBox(text + '\n\nSaved to ' + path, 'Mesh to Part - report')


def write_report(text):
    try:
        os.makedirs(REPORT_DIR, exist_ok=True)
        path = os.path.join(REPORT_DIR, time.strftime('report_%Y%m%d_%H%M%S.txt'))
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write(text)
        return path
    except Exception:
        return LOG_PATH


def on_validate(inputs):
    try:
        return (inputs.itemById('selMesh').selectionCount > 0
                and inputs.itemById('maxError').value > 0)
    except Exception:
        return False


# --------------------------------------------------------------------------- #
# dev bridge commands - how this add-in gets tested inside Fusion
# --------------------------------------------------------------------------- #

def dev_ping(app, job):
    des = adsk.fusion.Design.cast(app.activeProduct)
    bodies = []
    if des:
        found = []
        collect_mesh_bodies(des.rootComponent, found)
        for b in found:
            bb = b.boundingBox
            bodies.append({
                'name': b.name, 'tri': b.mesh.triangleCount, 'nodes': b.mesh.nodeCount,
                'closed': b.isClosed, 'oriented': b.isOriented,
                'size_mm': [round((bb.maxPoint.x - bb.minPoint.x) * CM, 3),
                            round((bb.maxPoint.y - bb.minPoint.y) * CM, 3),
                            round((bb.maxPoint.z - bb.minPoint.z) * CM, 3)],
            })
    return {
        'version': VERSION,
        'document': app.activeDocument.name if app.activeDocument else None,
        'design_type': (None if not des else
                        ('parametric' if is_parametric(des) else 'direct')),
        'brep_bodies': 0 if not des else des.rootComponent.bRepBodies.count,
        'mesh_bodies': bodies,
        'solid_bodies': [] if not des else collect_brep_info(des.rootComponent),
    }


def collect_brep_bodies(comp, out):
    """Every BRep body in the assembly, at any depth (see collect_mesh_bodies)."""
    for i in range(comp.bRepBodies.count):
        out.append(comp.bRepBodies.item(i))
    try:
        occs = comp.allOccurrences
    except Exception:
        return
    for i in range(occs.count):
        occ = occs.item(i)
        try:
            bodies = occ.component.bRepBodies
            for j in range(bodies.count):
                out.append(bodies.item(j).createForAssemblyContext(occ))
        except Exception:
            pass


def collect_brep_info(comp):
    found = []
    collect_brep_bodies(comp, found)
    out = []
    for b in found:
        try:
            bb = b.boundingBox
            out.append({
                'name': b.name,
                'component': b.parentComponent.name,
                'faces': b.faces.count,
                'solid': b.isSolid,
                'visible': b.isVisible,
                'volume_mm3': round(b.volume * 1000.0, 1),
                'size_mm': [round((bb.maxPoint.x - bb.minPoint.x) * CM, 2),
                            round((bb.maxPoint.y - bb.minPoint.y) * CM, 2),
                            round((bb.maxPoint.z - bb.minPoint.z) * CM, 2)],
            })
        except Exception:
            pass
    return out


def dev_probe(app, job):
    des = adsk.fusion.Design.cast(app.activeProduct)
    found = []
    collect_mesh_bodies(des.rootComponent, found)
    if not found:
        return {'error': 'no mesh body in the active document'}
    b = _pick_body(found, job)
    avail = probe_methods(b.parentComponent, b)
    return {'body': b.name,
            'methods': {k: {'ok': v[0], 'why': v[1]} for k, v in avail.items()}}


def dev_analyse(app, job):
    des = adsk.fusion.Design.cast(app.activeProduct)
    found = []
    collect_mesh_bodies(des.rootComponent, found)
    if not found:
        return {'error': 'no mesh body in the active document'}
    b = _pick_body(found, job)
    t0 = time.time()
    coords, idx = read_triangles(b)
    wc, wi = weld(coords, idx)
    sizes = sorted((loop_span_mm(wc, lp) for lp in boundary_loops(wi)), reverse=True)
    return {'body': b.name, 'tri': len(idx) // 3, 'welded_nodes': len(wc) // 3,
            'closed': b.isClosed, 'loops': [round(v, 4) for v in sizes[:40]],
            'n_loops': len(sizes), 'secs': round(time.time() - t0, 2)}


def test_sphere(n_lat=24, n_lon=36, radius_mm=15.0, hole_rows=2):
    """Synthetic UV sphere with an optional hole punched around the north pole.
    Returned in cm as (coords, idx). Kept out of the job handler so it can be unit
    tested outside Fusion."""
    r = float(radius_mm) * MM

    verts = [(0.0, 0.0, r)]
    for i in range(1, n_lat):
        th = math.pi * i / n_lat
        for j in range(n_lon):
            ph = 2.0 * math.pi * j / n_lon
            verts.append((r * math.sin(th) * math.cos(ph),
                          r * math.sin(th) * math.sin(ph),
                          r * math.cos(th)))
    verts.append((0.0, 0.0, -r))
    south = len(verts) - 1

    def ring(i, j):
        return 1 + (i - 1) * n_lon + (j % n_lon)

    # winding has to agree between the caps and the band, or every shared edge
    # reads as a boundary edge; the whole mesh is flipped afterwards if it came
    # out inside-out.
    tris = []
    for j in range(n_lon):
        tris.append((0, ring(1, j + 1), ring(1, j)))
    for i in range(1, n_lat - 1):
        for j in range(n_lon):
            a, b, c, d = ring(i, j), ring(i, j + 1), ring(i + 1, j + 1), ring(i + 1, j)
            tris.append((a, b, c))
            tris.append((a, c, d))
    for j in range(n_lon):
        tris.append((south, ring(n_lat - 1, j), ring(n_lat - 1, j + 1)))

    vol = 0.0
    for a, b, c in tris:
        pa, pb, pc = verts[a], verts[b], verts[c]
        vol += (pa[0] * (pb[1] * pc[2] - pb[2] * pc[1])
                - pa[1] * (pb[0] * pc[2] - pb[2] * pc[0])
                + pa[2] * (pb[0] * pc[1] - pb[1] * pc[0]))
    if vol < 0:
        tris = [(a, c, b) for a, b, c in tris]

    if hole_rows:                       # punch a hole around the north pole
        keep = []
        drop_rows = set(range(0, hole_rows + 1))
        for t in tris:
            rows = []
            for v in t:
                if v == 0:
                    rows.append(0)
                elif v == south:
                    rows.append(n_lat)
                else:
                    rows.append((v - 1) // n_lon + 1)
            if all(x in drop_rows for x in rows) and max(rows) <= hole_rows:
                continue
            keep.append(t)
        tris = keep

    coords, idx = [], []
    for v in verts:
        coords.extend(v)
    for t in tris:
        idx.extend(t)
    return coords, idx


def dev_make_test_doc(app, job):
    """A fresh document holding the synthetic sphere, so tests never touch the
    user's design."""
    coords, idx = test_sphere(int(job.get('lat', 24)), int(job.get('lon', 36)),
                              float(job.get('radius_mm', 15.0)),
                              int(job.get('hole_rows', 2)))
    doc = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
    des = adsk.fusion.Design.cast(app.activeProduct)
    # new documents follow the user's "capture design history" preference, which is
    # off here; Form features need a timeline, so force it on for scratch docs
    if not is_parametric(des):
        des.designType = adsk.fusion.DesignTypes.ParametricDesignType
    body, _base, _n = add_mesh_body(des.rootComponent, is_parametric(des),
                                    coords, idx, 'TestSphere')
    return {'document': doc.name, 'body': body.name if body else None,
            'tri': len(idx) // 3, 'nodes': len(coords) // 3,
            'closed': body.isClosed if body else None}


def dev_clone_to_new_doc(app, job):
    """Copy a mesh body from the active document into a brand new one, so tests
    never touch the user's design."""
    des = adsk.fusion.Design.cast(app.activeProduct)
    found = []
    collect_mesh_bodies(des.rootComponent, found)
    if not found:
        return {'error': 'no mesh body in the active document'}
    b = _pick_body(found, job)
    name = b.name
    coords, idx = read_triangles(b)
    doc = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
    des2 = adsk.fusion.Design.cast(app.activeProduct)
    if not is_parametric(des2):
        des2.designType = adsk.fusion.DesignTypes.ParametricDesignType
    body, _base, _n = add_mesh_body(des2.rootComponent, is_parametric(des2),
                                    coords, idx, name)
    return {'document': doc.name, 'body': body.name if body else None,
            'tri': len(idx) // 3}


def dev_convert(app, job):
    """Run the whole pipeline headlessly with an options dict."""
    des = adsk.fusion.Design.cast(app.activeProduct)
    found = []
    collect_mesh_bodies(des.rootComponent, found)
    # solids are convertible too now: the pipeline tessellates them at entry
    collect_brep_bodies(des.rootComponent, found)
    if not found:
        return {'error': 'no body in the active document'}
    b = _pick_body(found, job)
    o = dev_options(job.get('options', {}))
    t0 = time.time()
    rep = run_pipeline(app, b, o, None)
    rep['text'] = format_report(rep, o, time.time() - t0)
    rep['secs'] = round(time.time() - t0, 2)
    return rep


def dev_options(over):
    """Option dict for headless runs; same keys resolve_options produces."""
    o = {
        'max_error_mm': 0.1, 'escalate': True, 'max_attempts': 3,
        'fill_gaps': True, 'max_gap_mm': 1.0, 'fill_all': False,
        'repair_type': None, 'repair_label': 'None', 'rebuild_density': 128,
        'reduce': True, 'reduce_over': 150000, 'reduce_share': 0.4,
        'faceted_cap': 20000, 'fallback_faceted': True,
        'cage': 24, 'cage_auto': False, 'fairness': 0.01, 'ctrl_ratio': 1.0,
        'wrap_offset_mm': 0.0, 'wrap_close_mm': 0.0, 'wrap_smooth': 15,
        'aspect_limit': 3.0, 'aspect_iters': 80, 'aspect_gain': 2.0,
        'wrap_iters': 50, 'relax_iters': 70,
        'method_key': 'prismatic', 'resolution': 'Auto from max error',
        'face_count': 5000, 'preprocess_holes': True, 'face_groups': True,
        'feature_angle_rad': math.radians(25.0), 'min_group': 0.02,
        'op_type': adsk.fusion.MeshConvertOperationTypes
                   .ParametricFeatureMeshConvertOperationType,
        'verify': True, 'hide_source': True, 'show_report': False,
    }
    if 'feature_angle_deg' in over:
        o['feature_angle_rad'] = math.radians(float(over.pop('feature_angle_deg')))
    if 'repair' in over:
        label = over.pop('repair')
        for lab, en in REPAIR_MODES:
            if lab == label:
                o['repair_type'], o['repair_label'] = en, lab
    o.update(over)
    o['method_label'] = METHOD_BY_KEY[o['method_key']][0]
    return o


def dev_capture(app, job):
    vp = app.activeViewport
    cam = vp.camera
    cam.viewOrientation = getattr(adsk.core.ViewOrientations,
                                  job.get('view', 'IsoTopRightViewOrientation'))
    cam.isSmoothTransition = False
    vp.camera = cam
    adsk.doEvents()
    # fit() replaces the camera object, so the orientation has to be re-applied
    # after it or every shot comes out from whatever angle was already set
    vp.fit()
    adsk.doEvents()
    cam = vp.camera
    cam.viewOrientation = getattr(adsk.core.ViewOrientations,
                                  job.get('view', 'IsoTopRightViewOrientation'))
    cam.isSmoothTransition = False
    vp.camera = cam
    adsk.doEvents()
    vp.fit()
    adsk.doEvents()
    vp.refresh()
    path = job.get('path') or os.path.join(ADDIN_DIR, 'devjobs', 'capture.png')
    os.makedirs(os.path.dirname(path), exist_ok=True)
    ok = vp.saveAsImageFile(path, int(job.get('w', 900)), int(job.get('h', 700)))
    return {'ok': bool(ok), 'path': path}


def dev_close_doc(app, job):
    doc = app.activeDocument
    name = doc.name
    doc.close(False)
    return {'closed': name,
            'now': app.activeDocument.name if app.activeDocument else None}


def dev_panels(app, job):
    ui = app.userInterface
    ids = [ui.allToolbarPanels.item(i).id for i in range(ui.allToolbarPanels.count)]
    return {'mesh_panels': [p for p in ids if 'esh' in p],
            'solid_panels': [p for p in ids if p.startswith('Solid')][:20],
            'total': len(ids)}


CMD_ID = 'MeshToPart_Convert'

# Verified against this Fusion build: there is no plain "MeshModifyPanel"; the mesh
# tab panels are all ParaMesh*. Kept here rather than in the shell so the placement
# can be corrected without reloading the add-in.
PANEL_GROUPS = [
    ['ParaMeshModifyPanel', 'ParaMeshPreparePanel', 'ParaMeshCreatePanel',
     'MeshModifyPanel'],
    ['SolidCreatePanel'],
]


def place_button(ui, cmd_def):
    """Put the button in the first panel of each group that exists. Returns the ids."""
    placed = []
    for group in PANEL_GROUPS:
        for pid in group:
            panel = ui.allToolbarPanels.itemById(pid)
            if not panel:
                continue
            existing = panel.controls.itemById(cmd_def.id)
            if existing:
                existing.deleteMe()
            ctrl = panel.controls.addCommand(cmd_def)
            ctrl.isPromoted = True
            ctrl.isPromotedByDefault = True
            placed.append(pid)
            break
    return placed


def dev_rebuild_toolbar(app, job):
    ui = app.userInterface
    cd = ui.commandDefinitions.itemById(CMD_ID)
    if not cd:
        return {'error': 'command definition %s not found' % CMD_ID}
    return {'placed': place_button(ui, cd)}


def dev_purge_work(app, job):
    """Delete leftover [work] and throw-away bodies from failed runs."""
    des = adsk.fusion.Design.cast(app.activeProduct)
    comp = des.rootComponent
    removed = []
    # collect_mesh_bodies, not comp.meshBodies: the work meshes live wherever
    # their source did, which for an imported STL assembly is inside a
    # component. Scanning only the root found nothing and quietly reported
    # success while the leftovers piled up.
    found = []
    collect_mesh_bodies(comp, found)
    for b in reversed(found):
        try:
            name_probe = b.name
        except Exception:
            continue
        if '[work' not in name_probe and not name_probe.startswith(TMP_NAME):
            continue
        name = name_probe
        if not job.get('confirm'):
            removed.append(name)          # dry run: report, change nothing
            continue
        try:
            b.deleteMe()                  # the BODY only, never its BaseFeature
            removed.append(name)
        except Exception:
            log('could not delete ' + name)
            log(traceback.format_exc())

    rest = []
    collect_mesh_bodies(comp, rest)
    return {'dry_run': not bool(job.get('confirm')),
            'removed' if job.get('confirm') else 'would_remove': removed,
            'remaining': len(rest)}


def dev_open_dialog(app, job):
    """Open the real command dialog. on_created writes dialog_dump.json as it
    builds, which is how the UI path gets tested instead of eyeballed."""
    try:
        if os.path.exists(DIALOG_DUMP):
            os.remove(DIALOG_DUMP)
    except Exception:
        pass
    cd = app.userInterface.commandDefinitions.itemById(CMD_ID)
    if not cd:
        return {'error': 'command definition %s not found' % CMD_ID}
    return {'executed': bool(cd.execute())}


def dev_read_dialog(app, job):
    """Read the dump the dialog wrote, then close the dialog again."""
    ui = app.userInterface
    active = ui.activeCommand
    closed = False
    try:
        closed = bool(ui.terminateActiveCommand())
    except Exception:
        log(traceback.format_exc())
    dump = None
    try:
        with open(DIALOG_DUMP, encoding='utf-8') as fh:
            dump = json.load(fh)
    except Exception:
        pass
    return {'active_command': active, 'closed': closed, 'dump': dump}


def dev_press_ok(app, job):
    """Press OK on the open dialog - the full user path, execute handler included.
    Only ever aimed at a scratch document."""
    cmd = getattr(_STATE, 'live_cmd', None)
    if cmd is None:
        return {'error': 'no dialog has been opened in this session'}
    ok = cmd.doExecute(True)
    return {'executed': bool(ok), 'document': app.activeDocument.name}


def dev_set_inputs(app, job):
    """Set dialog inputs before pressing OK. values is {input_id: value}."""
    cmd = getattr(_STATE, 'live_cmd', None)
    if cmd is None:
        return {'error': 'no dialog open'}
    inputs = cmd.commandInputs
    done = {}
    for k, v in (job.get('values') or {}).items():
        ci = inputs.itemById(k)
        if ci is None:
            done[k] = 'no such input'
            continue
        try:
            if hasattr(ci, 'listItems') and isinstance(v, str):
                for n in range(ci.listItems.count):
                    if ci.listItems.item(n).name.startswith(v):
                        ci.listItems.item(n).isSelected = True
                        done[k] = ci.selectedItem.name
                        break
                else:
                    done[k] = 'no item starting with %r' % v
            elif isinstance(v, str) and hasattr(ci, 'expression'):
                ci.expression = v          # e.g. "0.1 mm"; avoids guessing cm
                done[k] = '%s -> %s' % (ci.expression, ci.value)
            else:
                ci.value = v
                done[k] = ci.value
        except Exception:
            done[k] = short_error(traceback.format_exc())
    try:
        refresh_visibility(inputs)
    except Exception:
        log(traceback.format_exc())
    return {'set': done}


def dev_tspline_probe(app, job):
    """Can a Form body be made at all on this licence, and what does the TSM parser
    actually want? Everything about a custom Form path depends on these answers."""
    out = {}
    ui = app.userInterface
    des = adsk.fusion.Design.cast(app.activeProduct)
    comp = des.rootComponent
    out['forms_before'] = comp.features.formFeatures.count

    try:
        ff = comp.features.formFeatures.add()
        out['formFeatures.add'] = 'ok'
        out['tsplineBodies_count'] = ff.tSplineBodies.count
    except Exception:
        out['formFeatures.add'] = short_error(traceback.format_exc())
        ff = None

    if ff:
        probes = job.get('tsm') or ['', 'tsm 1', 'TSM']
        out['addByTSMDescription'] = {}
        for t in probes:
            try:
                b = ff.tSplineBodies.addByTSMDescription(t)
                out['addByTSMDescription'][t[:40]] = 'ok: %s' % (b.name if b else None)
            except Exception:
                out['addByTSMDescription'][t[:40]] = short_error(traceback.format_exc())
        try:
            ff.deleteMe()
            out['form_deleted'] = True
        except Exception:
            out['form_deleted'] = short_error(traceback.format_exc())

    names = []
    for i in range(ui.commandDefinitions.count):
        cd = ui.commandDefinitions.item(i)
        try:
            cid = cd.id
        except Exception:
            continue
        low = cid.lower()
        if 'tspline' in low or 'sculpt' in low or 'form' in low or 'quadball' in low:
            names.append(cid)
    out['commands'] = sorted(names)[:120]
    return out


def dev_tsm_try(app, job):
    """Feed candidate TSM strings to the parser inside an ACTIVE form feature and
    report what it says. The error text is the only documentation available."""
    des = adsk.fusion.Design.cast(app.activeProduct)
    comp = des.rootComponent
    out = {}
    for name, text in (job.get('cases') or {}).items():
        ff = None
        try:
            ff = comp.features.formFeatures.add()
            ff.startEdit()
            b = ff.tSplineBodies.addByTSMDescription(text)
            out[name] = {'ok': True, 'body': b.name if b else None,
                         'count': ff.tSplineBodies.count}
            try:
                out[name]['roundtrip'] = b.getTSMDescription()[:2000] if b else None
            except Exception:
                pass
        except Exception:
            out[name] = {'ok': False, 'why': short_error(traceback.format_exc())}
        finally:
            try:
                if ff:
                    ff.finishEdit()
                    ff.deleteMe()
            except Exception:
                log(traceback.format_exc())
    return out


def dev_arm_primitive(app, job):
    """Hook Fusion's own T-Spline primitive command so it can be driven from here,
    which is how the real TSM format gets learned."""
    ui = app.userInterface
    cid = job.get('command', 'TSplinePrimitiveQuadBallCommand')
    cd = ui.commandDefinitions.itemById(cid)
    if not cd:
        return {'error': 'no command %s' % cid}

    class _Grab(adsk.core.CommandCreatedEventHandler):
        def notify(self, args):
            try:
                _STATE.prim_cmd = args.command
            except Exception:
                log(traceback.format_exc())

    h = _Grab()
    _STATE.prim_handler = h          # must outlive this function
    cd.commandCreated.add(h)
    _STATE.prim_def = cd
    _STATE.prim_cmd = None

    out = {'command': cid}
    if job.get('enter_form'):
        # a T-Spline primitive only runs inside an active form feature
        des = adsk.fusion.Design.cast(app.activeProduct)
        ff = des.rootComponent.features.formFeatures.add()
        ff.startEdit()
        _STATE.form_feature = ff
        out['form_started'] = True
    out['executed'] = bool(cd.execute())
    out['activeCommand'] = ui.activeCommand
    return out


def dev_fire_primitive(app, job):
    """Press OK on the armed primitive command, then dump whatever T-Spline it made."""
    cmd = getattr(_STATE, 'prim_cmd', None)
    out = {'had_command': cmd is not None}
    if cmd is not None:
        try:
            out['inputs'] = [cmd.commandInputs.item(i).id
                             for i in range(cmd.commandInputs.count)]
        except Exception:
            out['inputs'] = short_error(traceback.format_exc())
        try:
            out['doExecute'] = bool(cmd.doExecute(True))
        except Exception:
            out['doExecute'] = short_error(traceback.format_exc())
    out['activeCommand'] = app.userInterface.activeCommand
    try:
        app.userInterface.terminateActiveCommand()
    except Exception:
        pass
    ff = getattr(_STATE, 'form_feature', None)
    if ff is not None:
        try:
            ff.finishEdit()
            out['form_finished'] = True
        except Exception:
            out['form_finished'] = short_error(traceback.format_exc())
        _STATE.form_feature = None

    des = adsk.fusion.Design.cast(app.activeProduct)
    comp = des.rootComponent
    out['form_features'] = comp.features.formFeatures.count
    dumps = []
    for i in range(comp.features.formFeatures.count):
        ff = comp.features.formFeatures.item(i)
        for k in range(ff.tSplineBodies.count):
            tb = ff.tSplineBodies.item(k)
            try:
                dumps.append({'name': tb.name, 'tsm': tb.getTSMDescription()})
            except Exception:
                dumps.append({'name': tb.name,
                              'error': short_error(traceback.format_exc())})
    out['bodies'] = dumps
    if dumps and job.get('save'):
        path = os.path.join(ADDIN_DIR, 'devjobs', 'sample.tsm')
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write(dumps[0].get('tsm') or '')
        out['saved'] = path
    return out


NO_WINDOW = 0x08000000 if os.name == 'nt' else 0   # or a console flashes per call

PYTHON_CANDIDATES = [
    os.path.expandvars(r'%LOCALAPPDATA%\Programs\Python\Python312\python.exe'),
    os.path.expandvars(r'%LOCALAPPDATA%\Programs\Python\Python311\python.exe'),
    os.path.expandvars(r'%LOCALAPPDATA%\Programs\Python\Python310\python.exe'),
    r'C:\Python312\python.exe',
    r'C:\Python311\python.exe',
    r'C:\Python310\python.exe',
    '/usr/bin/python3',
    '/usr/local/bin/python3',
]


def find_python():
    """Locate an interpreter that has numpy, scipy and trimesh.

    Fusion's own Python has none of them, so the fitting runs outside. Order:
    the MESHTOPART_PYTHON environment variable, a python_path.txt beside the
    add-in, the usual install locations, then the py launcher.
    """
    import subprocess

    def works(exe, args=None):
        try:
            cmd = ([exe] + (args or []) +
                   ['-c', 'import numpy, scipy, trimesh'])
            r = subprocess.run(cmd, capture_output=True, timeout=60,
                               creationflags=NO_WINDOW)
            return r.returncode == 0
        except Exception:
            return False

    env = os.environ.get('MESHTOPART_PYTHON')
    if env and works(env):
        return [env]

    cfg = os.path.join(ADDIN_DIR, 'python_path.txt')
    if os.path.exists(cfg):
        try:
            with open(cfg, encoding='utf-8') as fh:
                p = fh.read().strip()
            if p and works(p):
                return [p]
        except Exception:
            pass

    for p in PYTHON_CANDIDATES:
        if os.path.exists(p) and works(p):
            return [p]

    for args in (['-3.12'], ['-3.11'], ['-3.10'], ['-3']):
        if works('py', args):
            return ['py'] + args
    return None


_PYTHON = None


def python_cmd():
    global _PYTHON
    if _PYTHON is None:
        _PYTHON = find_python() or []
    return list(_PYTHON)


def dev_subprocess_test(app, job):
    """Fusion's bundled Python has no numpy. The real fitting has to run in py -3.10,
    so the whole organic path depends on being able to launch it from in here."""
    import subprocess
    exe = python_cmd()
    out = {'interpreter': exe}
    if not exe:
        return out
    try:
        r = subprocess.run(exe + ['-c',
                            'import numpy, scipy, trimesh, sys;'
                            'print(numpy.__version__, scipy.__version__,'
                            ' trimesh.__version__, sys.version.split()[0])'],
                           capture_output=True, text=True, timeout=120,
                           creationflags=NO_WINDOW)
        out['returncode'] = r.returncode
        out['stdout'] = r.stdout.strip()
        out['stderr'] = r.stderr.strip()[:500]
    except Exception:
        out['error'] = short_error(traceback.format_exc())
    return out


def dev_export_mesh(app, job):
    """Dump a mesh body's triangles as JSON so the py -3.10 engine can work on it."""
    des = adsk.fusion.Design.cast(app.activeProduct)
    found = []
    collect_mesh_bodies(des.rootComponent, found)
    if not found:
        return {'error': 'no mesh body in the active document'}
    b = _pick_body(found, job)
    coords, idx = read_triangles(b)
    path = job.get('path') or os.path.join(ADDIN_DIR, 'devjobs', 'mesh.json')
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump({'coords': coords, 'idx': idx, 'name': b.name}, fh)
    return {'path': path, 'body': b.name, 'tri': len(idx) // 3}


# --------------------------------------------------------------------------- #
# custom organic conversion - our own replacement for the gated Organic method
#
# The heavy fitting runs in py -3.10 (numpy/scipy/trimesh); this side turns the
# six fitted NURBS patches into a real B-Rep solid through BRepBodyDefinition.
# --------------------------------------------------------------------------- #

ORGANIC_ENGINE = os.path.join(ADDIN_DIR, 'mtp_organic.py')
WRAP_ENGINE = os.path.join(ADDIN_DIR, 'mtp_wrap.py')


def _pt(v):
    return adsk.core.Point3D.create(v[0], v[1], v[2])


def run_engine(coords, idx, params, timeout=1800):
    """Hand the mesh to the fitting engine and get patches + topology back."""
    import subprocess
    work = os.path.join(ADDIN_DIR, 'devjobs')
    os.makedirs(work, exist_ok=True)
    jp = os.path.join(work, 'organic_job.json')
    op = os.path.join(work, 'organic_out.json')
    job = dict(params)
    job['coords'] = coords
    job['idx'] = idx
    with open(jp, 'w', encoding='utf-8') as fh:
        json.dump(job, fh)
    if os.path.exists(op):
        os.remove(op)
    exe = python_cmd()
    if not exe:
        raise RuntimeError(
            'No Python with numpy, scipy and trimesh was found. Set the '
            'MESHTOPART_PYTHON environment variable, or put the full path to the '
            'interpreter in python_path.txt next to the add-in.')
    r = subprocess.run(exe + [ORGANIC_ENGINE, jp, op], capture_output=True,
                       text=True, timeout=timeout, creationflags=NO_WINDOW)
    if not os.path.exists(op):
        raise RuntimeError('fitting engine produced nothing (rc=%s)\n%s\n%s'
                           % (r.returncode, r.stdout[-800:], r.stderr[-800:]))
    with open(op, encoding='utf-8') as fh:
        res = json.load(fh)
    if not res.get('ok'):
        raise RuntimeError('fitting engine failed:\n%s'
                           % res.get('error', '')[-1200:])
    return res


def surface_from_patch(p, swap):
    """NurbsSurface from the engine's control net.

    The engine emits the net as [u][v] with u major. Fusion's expected ordering
    for the flat control point array is not documented, so it is established once
    by round-trip and passed in here rather than assumed.
    """
    cu, cv = int(p['countU']), int(p['countV'])
    pts = p['controlPoints']
    if swap:
        order = [pts[u * cv + v] for v in range(cv) for u in range(cu)]
    else:
        order = pts
    return adsk.core.NurbsSurface.create(
        3, 3, cu, cv, [_pt(q) for q in order],
        p['knotsU'], p['knotsV'], [],          # non-rational: no weights
        adsk.core.NurbsSurfaceProperties.OpenNurbsSurface,
        adsk.core.NurbsSurfaceProperties.OpenNurbsSurface)


def detect_cp_order():
    """Work out whether Fusion wants the control net u-major or v-major.

    Build a deliberately asymmetric 4x5 net, evaluate the surface near one
    corner, and see which interpretation the answer agrees with. Guessing this
    silently transposes every patch.
    """
    cu, cv = 4, 5
    net = [[float(u), float(v), float(u * v) * 0.1]
           for u in range(cu) for v in range(cv)]        # [u][v], u major
    ku = [0, 0, 0, 0, 1, 1, 1, 1]
    kv = [0, 0, 0, 0, 0.5, 1, 1, 1, 1]
    surf = adsk.core.NurbsSurface.create(
        3, 3, cu, cv, [_pt(q) for q in net], ku, kv, [],
        adsk.core.NurbsSurfaceProperties.OpenNurbsSurface,
        adsk.core.NurbsSurfaceProperties.OpenNurbsSurface)
    if surf is None:
        return False, 'NurbsSurface.create returned null'
    ev = surf.evaluator
    ok, p0 = ev.getPointAtParameter(adsk.core.Point2D.create(0.0, 0.0))
    if not ok:
        return False, 'evaluator failed'
    # u-major: corner (0,0) is net[0] = (0,0,0); v-major it is also (0,0,0),
    # so probe the far corner in u instead
    ok, p1 = ev.getPointAtParameter(adsk.core.Point2D.create(1.0, 0.0))
    umajor = abs(p1.x - 3.0) < 1e-6 and abs(p1.y - 0.0) < 1e-6
    vmajor = abs(p1.x - 0.0) < 1e-6 and abs(p1.y - 4.0) < 1e-6
    return (not umajor), 'u_major=%s v_major=%s corner=(%.3f, %.3f, %.3f)' % (
        umajor, vmajor, p1.x, p1.y, p1.z)


def build_brep(res, swap, log_lines):
    """Assemble the six patches into one solid: 8 vertices, 12 edges, 6 faces."""
    topo = res['topology']
    bd = adsk.fusion.BRepBodyDefinition.create()
    bd.doFullHealing = True
    lump = bd.lumpDefinitions.add()
    shell = lump.shellDefinitions.add()

    vdefs = {}
    for vid, pos in topo['vertices'].items():
        vdefs[vid] = bd.createVertexDefinition(_pt(pos))

    edefs = {}
    for key, e in topo['edges'].items():
        curve = adsk.core.NurbsCurve3D.createNonRational(
            [_pt(q) for q in e['controlPoints']], 3, e['knots'], False)
        if curve is None:
            raise RuntimeError('edge %s: NurbsCurve3D.createNonRational failed' % key)
        edefs[key] = bd.createEdgeDefinitionByCurve(
            vdefs[str(e['start'])], vdefs[str(e['end'])], curve)

    for f, loop_spec in enumerate(topo['faces']):
        surf = surface_from_patch(res['patches'][f], swap)
        if surf is None:
            raise RuntimeError('patch %d: NurbsSurface.create failed' % f)
        fd = shell.faceDefinitions.add(surf, False)
        loop = fd.loopDefinitions.add()
        for item in loop_spec:
            loop.bRepCoEdgeDefinitions.add(edefs[item['edge']], item['reversed'])

    body = bd.createBody()
    log_lines.extend(list(bd.outcomeInfo or []))
    return body


def read_brep_triangles(body):
    """(coords_cm, tri_idx) from a SOLID body, via its display tessellation.

    displayMeshes.bestMesh, not a TriangleMeshCalculator: the calculator can
    take Fusion down on a large organic body, and the display mesh is already
    as fine as the viewport is showing.
    """
    tm = body.meshManager.displayMeshes.bestMesh
    return list(tm.nodeCoordinatesAsDouble), list(tm.nodeIndices)


def run_wrap(coords, idx, params, timeout=3600):
    """Hand a mesh to the offset engine and get the offset mesh back."""
    import subprocess
    work = os.path.join(ADDIN_DIR, 'devjobs')
    os.makedirs(work, exist_ok=True)
    jp = os.path.join(work, 'wrap_job.json')
    op = os.path.join(work, 'wrap_out.json')
    job = dict(params)
    job['coords'] = coords
    job['idx'] = idx
    with open(jp, 'w', encoding='utf-8') as fh:
        json.dump(job, fh)
    if os.path.exists(op):
        os.remove(op)
    exe = python_cmd()
    if not exe:
        raise RuntimeError(
            'No Python with numpy, scipy, trimesh and pyvista was found. Set '
            'MESHTOPART_PYTHON, or put the interpreter path in python_path.txt '
            'next to the add-in.')
    r = subprocess.run(exe + [WRAP_ENGINE, jp, op], capture_output=True,
                       text=True, timeout=timeout, creationflags=NO_WINDOW)
    if not os.path.exists(op):
        raise RuntimeError('wrap engine produced nothing (rc=%s)\n%s\n%s'
                           % (r.returncode, r.stdout[-800:], r.stderr[-800:]))
    with open(op, encoding='utf-8') as fh:
        res = json.load(fh)
    if not res.get('ok'):
        raise RuntimeError('wrap engine failed:\n%s' % res.get('error', '')[-1200:])
    return res


def wrap_body(app, src_body, params, log_lines):
    """Offset any body -- mesh or solid -- and drop the result in as a mesh.

    Returns (mesh_body, engine result). Fusion's Offset Face and Shell are per
    face and give up when an outward offset self-intersects in a concavity;
    this goes through a signed distance field, where that case simply cannot
    arise.
    """
    design = adsk.fusion.Design.cast(app.activeProduct)
    parametric = is_parametric(design)
    comp = src_body.parentComponent

    if adsk.fusion.BRepBody.cast(src_body):
        coords, idx = read_brep_triangles(src_body)
    else:
        coords, idx = read_triangles(src_body)

    res = run_wrap(coords, idx, params)
    log_lines.extend(res.get('log', []))

    name = '%s [wrap %+.2f]' % (src_body.name, params['offset_mm'])
    body, base, actual = add_mesh_body(comp, parametric, res['coords'],
                                       res['idx'], name)
    res['stats']['body'] = actual
    return body, res, base


def organic_body(app, src_body, params, log_lines):
    """Fit and build. Returns (BRepBody, engine result)."""
    design = adsk.fusion.Design.cast(app.activeProduct)
    parametric = is_parametric(design)
    comp = src_body.parentComponent

    coords, idx = read_triangles(src_body)
    res = run_engine(coords, idx, params)
    log_lines.extend(res.get('log', []))

    swap, note = detect_cp_order()
    log_lines.append('control point order: %s' % note)

    temp = build_brep(res, swap, log_lines)
    if temp is None:
        raise RuntimeError('BRepBodyDefinition.createBody failed:\n  %s'
                           % '\n  '.join(log_lines[-8:]))

    base = None
    if parametric:
        base = comp.features.baseFeatures.add()
        base.startEdit()
    try:
        body = comp.bRepBodies.add(temp, base)
    finally:
        if base:
            base.finishEdit()
    return body, res, base


def dev_organic(app, job):
    """End-to-end custom organic conversion, measured against the source mesh."""
    des = adsk.fusion.Design.cast(app.activeProduct)
    found = []
    collect_mesh_bodies(des.rootComponent, found)
    if not found:
        return {'error': 'no mesh body in the active document'}
    src = _pick_body(found, job)
    comp = src.parentComponent
    log_lines = []
    t0 = time.time()

    params = {'cage': 24, 'wrap_iters': 50, 'relax_iters': 60,
              'ctrl': 25, 'fair': 1e-2}
    params.update(job.get('params') or {})

    body, res, base = organic_body(app, src, params, log_lines)
    body.name = '%s [organic]' % src.name

    out = {
        'body': body.name,
        'faces': body.faces.count,
        'edges': body.edges.count,
        'vertices': body.vertices.count,
        'is_solid': body.isSolid,
        'volume_mm3': body.volume * 1000.0 if body.isSolid else None,
        'src_volume_mm3': src.volume * 1000.0 if src.isClosed else None,
        'engine': res['stats'],
        'log': log_lines[-24:],
        'secs': round(time.time() - t0, 1),
    }
    if job.get('verify', True):
        dev = measure_deviation(comp, is_parametric(des), src, body,
                                float(job.get('max_error_mm', 1.0)))
        out['measured'] = None if not dev else {
            'part_to_mesh': {k: round(v, 4) for k, v in dev['out'].items()},
            'mesh_to_part': {k: round(v, 4) for k, v in dev['back'].items()},
        }
    if job.get('cleanup'):
        try:
            (base or body).deleteMe()
            out['cleaned'] = True
        except Exception:
            out['cleaned'] = short_error(traceback.format_exc())
    return out


def dev_cp_order(app, job):
    swap, note = detect_cp_order()
    return {'swap': swap, 'note': note}


def dev_read_tsm(app, job):
    """Dump the TSM description of any T-Spline body in the open documents.

    The format is undocumented and the parser gives no hints, so one hand-made
    Form primitive is all it takes to learn it - this is what reads it.
    """
    out = {'documents': [], 'bodies': []}
    for i in range(app.documents.count):
        doc = app.documents.item(i)
        out['documents'].append(doc.name)
    des = adsk.fusion.Design.cast(app.activeProduct)
    if not des:
        return out
    comp = des.rootComponent
    for i in range(comp.features.formFeatures.count):
        ff = comp.features.formFeatures.item(i)
        for k in range(ff.tSplineBodies.count):
            tb = ff.tSplineBodies.item(k)
            rec = {'name': tb.name, 'form_index': i}
            try:
                tsm = tb.getTSMDescription()
                rec['length'] = len(tsm)
                rec['head'] = tsm[:1500]
                path = os.path.join(ADDIN_DIR, 'devjobs', 'sample.tsm')
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with open(path, 'w', encoding='utf-8', errors='replace') as fh:
                    fh.write(tsm)
                rec['saved'] = path
            except Exception:
                rec['error'] = short_error(traceback.format_exc())
            out['bodies'].append(rec)
    out['form_features'] = comp.features.formFeatures.count
    return out


def dev_list_brep(app, job):
    """Name every BRep body, and optionally remove ones this tool created."""
    des = adsk.fusion.Design.cast(app.activeProduct)
    comp = des.rootComponent
    rows = []
    for i in range(comp.bRepBodies.count):
        b = comp.bRepBodies.item(i)
        rows.append({'name': b.name, 'solid': b.isSolid, 'faces': b.faces.count,
                     'volume_mm3': round(b.volume * 1000.0, 1) if b.isSolid else None})
    removed = []
    if job.get('remove_ours'):
        for i in range(comp.bRepBodies.count - 1, -1, -1):
            b = comp.bRepBodies.item(i)
            names = job.get('names')
            if names is not None:
                if b.name not in names:
                    continue
            elif '[organic]' not in b.name and '[part]' not in b.name:
                continue
            name = b.name
            try:
                b.deleteMe()          # BRepBody has no baseOrFormFeature
                removed.append(name)
            except Exception:
                log(traceback.format_exc())
    return {'bodies': rows, 'removed': removed,
            'remaining': comp.bRepBodies.count}


def dev_import_tsm(app, job):
    """Create a Form body from a TSM description. Reports the body it produced so
    the unit scale and orientation can be measured rather than assumed."""
    if job.get('new_doc'):
        app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
    des = adsk.fusion.Design.cast(app.activeProduct)
    if not is_parametric(des):
        des.designType = adsk.fusion.DesignTypes.ParametricDesignType
    comp = des.rootComponent
    text = job.get('text')
    if not text:
        with open(job['path'], encoding='utf-8') as fh:
            text = fh.read()
    if not is_parametric(des):
        des.designType = adsk.fusion.DesignTypes.ParametricDesignType
    ff = comp.features.formFeatures.add()
    ff.startEdit()
    try:
        tb = ff.tSplineBodies.addByTSMDescription(text)
    finally:
        ff.finishEdit()
    if tb is None:
        return {'error': 'addByTSMDescription returned null'}
    out = {'tspline_body': tb.name, 'form_features': comp.features.formFeatures.count,
           'tspline_count': ff.tSplineBodies.count}
    for i in range(comp.bRepBodies.count):
        b = comp.bRepBodies.item(i)
        if b.name == tb.name:
            bb = b.boundingBox
            out['brep'] = {
                'name': b.name, 'faces': b.faces.count, 'solid': b.isSolid,
                'volume_mm3': round(b.volume * 1000.0, 4) if b.isSolid else None,
                'size_mm': [round((bb.maxPoint.x - bb.minPoint.x) * CM, 4),
                            round((bb.maxPoint.y - bb.minPoint.y) * CM, 4),
                            round((bb.maxPoint.z - bb.minPoint.z) * CM, 4)],
                'centre_mm': [round((bb.maxPoint.x + bb.minPoint.x) * CM / 2, 4),
                              round((bb.maxPoint.y + bb.minPoint.y) * CM / 2, 4),
                              round((bb.maxPoint.z + bb.minPoint.z) * CM / 2, 4)],
            }
    if job.get('cleanup'):
        try:
            ff.deleteMe()
            out['cleaned'] = True
        except Exception:
            out['cleaned'] = short_error(traceback.format_exc())
    return out


def organic_form_body(app, src_body, params, log_lines):
    """Fit a T-Spline control cage and hand it to Fusion as a real Form body.

    The cage is solved so its Catmull-Clark LIMIT surface lands on the mesh, then
    serialised as TSM. Returns (BRepBody, engine result, FormFeature).
    """
    comp = src_body.parentComponent
    coords, idx = read_triangles(src_body)
    p = dict(params)
    p['output'] = 'tsm'
    res = run_engine(coords, idx, p)
    log_lines.extend(res.get('log', []))

    des = adsk.fusion.Design.cast(app.activeProduct)
    try:
        # a form feature cannot be created with the marker parked inside the base
        # feature that holds the working mesh
        des.timeline.moveToEnd()
    except Exception:
        log_lines.append(short_error(traceback.format_exc()))

    # A cage Fusion rejects (self-intersection) leaves the form feature half built,
    # and every later formFeatures.add() then fails with InternalValidationError.
    # Tearing it down here is what lets the ladder try a coarser cage.
    ff = comp.features.formFeatures.add()
    try:
        ff.startEdit()
        tb = ff.tSplineBodies.addByTSMDescription(res['tsm'])
        ff.finishEdit()
    except Exception:
        for step in (ff.finishEdit, ff.deleteMe):
            try:
                step()
            except Exception:
                pass
        raise
    if tb is None:
        try:
            ff.deleteMe()
        except Exception:
            pass
        raise RuntimeError('addByTSMDescription rejected the generated cage')

    body = comp.bRepBodies.itemByName(tb.name)
    if body is None:
        body = newest_brep(comp)
    return body, res, ff


def dev_visibility(app, job):
    """Show or hide bodies by name fragment, so documentation shots can be set up."""
    des = adsk.fusion.Design.cast(app.activeProduct)
    comp = des.rootComponent
    changed = []
    want = job.get('show')
    for coll in (comp.meshBodies, comp.bRepBodies):
        for i in range(coll.count):
            b = coll.item(i)
            on = (want is None) or (want in b.name)
            hide = job.get('hide')
            if hide and hide in b.name:
                on = False
            try:
                b.isLightBulbOn = on
                changed.append('%s=%s' % (b.name, on))
            except Exception:
                pass
    return {'bodies': changed}


def dev_purge_forms(app, job):
    """Delete form features from the ACTIVE document. Used to clean up test bodies;
    it never touches another document."""
    des = adsk.fusion.Design.cast(app.activeProduct)
    comp = des.rootComponent
    doc = app.activeDocument.name
    if job.get('only_document') and job['only_document'] != doc:
        return {'skipped': doc}
    removed = []
    for i in range(comp.features.formFeatures.count - 1, -1, -1):
        ff = comp.features.formFeatures.item(i)
        names = [ff.tSplineBodies.item(k).name
                 for k in range(ff.tSplineBodies.count)]
        try:
            ff.deleteMe()
            removed.append(names or ['<empty>'])
        except Exception:
            log(traceback.format_exc())
    return {'document': doc, 'removed': removed,
            'form_features_left': comp.features.formFeatures.count,
            'brep_left': comp.bRepBodies.count}


def _pick_body(found, job):
    want = job.get('body')
    if want:
        for b in found:
            if b.name == want:
                return b
        # falling through to found[0] here meant a mistyped or solid-only name
        # silently ran against a completely different body, and the failure that
        # followed pointed at the wrong geometry
        raise RuntimeError('body %r not found. Available: %s'
                           % (want, ', '.join(b.name for b in found) or 'none'))
    return found[0]


def dev_wrap(app, job):
    """Offset a body by a distance, optionally refitting it back to a solid."""
    des = adsk.fusion.Design.cast(app.activeProduct)
    want = job.get('body')
    found = []
    collect_brep_bodies(des.rootComponent, found)
    collect_mesh_bodies(des.rootComponent, found)
    src = None
    for b in found:
        if b.name == want:
            src = b
            break
    if src is None:
        return {'error': 'body %r not found' % want,
                'available': [b.name for b in found]}

    o = job.get('options', {})
    params = {'offset_mm': float(o.get('offset_mm', 2.0)),
              'pitch_mm': float(o.get('pitch_mm', 0.0)),
              'close_mm': float(o.get('close_mm', 0.0)),
              'smooth_iters': int(o.get('smooth_iters', 15)),
              'budget': float(o.get('budget', 15.0e6))}
    logs = []
    t0 = time.time()
    mesh_body, res, _base = wrap_body(app, src, params, logs)
    out = {'source': src.name, 'wrap': res['stats']}

    if o.get('solid'):
        cage = int(o.get('cage', 48))
        m = cage + 1
        fit = {'cage': cage,
               'wrap_iters': 50,
               'relax_iters': int(o.get('relax_iters', 200)),
               'ctrl': max(6, min(m, int(round(float(o.get('ctrl_ratio', 0.5)) * m)))),
               'fair': float(o.get('fairness', 0.01))}
        body, eres, _b2 = organic_body(app, mesh_body, fit, logs)
        try:
            body.name = '%s [wrap %+.2f solid]' % (src.name, params['offset_mm'])
        except Exception:
            pass
        out['solid'] = body.name
        out['fit'] = eres['stats']

    out['secs'] = round(time.time() - t0, 2)
    out['log'] = logs[-40:]
    return out


DEV_COMMANDS = {
    'ping': dev_ping,
    'probe': dev_probe,
    'panels': dev_panels,
    'analyse': dev_analyse,
    'make_test_doc': dev_make_test_doc,
    'clone_to_new_doc': dev_clone_to_new_doc,
    'convert': dev_convert,
    'capture': dev_capture,
    'close_doc': dev_close_doc,
    'rebuild_toolbar': dev_rebuild_toolbar,
    'purge_work': dev_purge_work,
    'wrap': dev_wrap,
    'open_dialog': dev_open_dialog,
    'read_dialog': dev_read_dialog,
    'set_inputs': dev_set_inputs,
    'press_ok': dev_press_ok,
    'tspline_probe': dev_tspline_probe,
    'tsm_try': dev_tsm_try,
    'subprocess_test': dev_subprocess_test,
    'export_mesh': dev_export_mesh,
    'organic': dev_organic,
    'cp_order': dev_cp_order,
    'read_tsm': dev_read_tsm,
    'list_brep': dev_list_brep,
    'import_tsm': dev_import_tsm,
    'visibility': dev_visibility,
    'purge_forms': dev_purge_forms,
    'arm_primitive': dev_arm_primitive,
    'fire_primitive': dev_fire_primitive,
}
