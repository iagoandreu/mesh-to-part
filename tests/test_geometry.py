"""Offline tests for the parts of the add-in that need no Fusion.

    py -3.10 tests/test_geometry.py

adsk is stubbed (see tests/stub), so mtp_core imports outside Fusion. Everything
checked here has an answer known in advance.
"""

import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'stub'))
sys.path.insert(0, os.path.dirname(HERE))

import mtp_core as M                                        # noqa: E402

FAILED = []


def check(name, cond, detail=''):
    print('  %s %s%s' % ('ok  ' if cond else 'FAIL', name,
                         '' if cond else '   %s' % (detail,)))
    if not cond:
        FAILED.append(name)


# --------------------------------------------------------------------------- #
# welding and boundary detection, on a box whose answers are known by hand
# --------------------------------------------------------------------------- #
S = 1.0
V = [(0, 0, 0), (S, 0, 0), (S, S, 0), (0, S, 0),
     (0, 0, S), (S, 0, S), (S, S, S), (0, S, S)]
FACES = [(0, 3, 2), (0, 2, 1), (4, 5, 6), (4, 6, 7),
         (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
         (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7)]


def build(tris):
    coords, idx = [], []
    for t in tris:
        for v in t:
            idx.append(len(coords) // 3)
            coords.extend(V[v])
    return coords, idx


def normal(a, b, c):
    u = [b[k] - a[k] for k in range(3)]
    v = [c[k] - a[k] for k in range(3)]
    return (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0])


print('welding and boundaries')
c, i = build(FACES)
wc, wi = M.weld(c, i)
check('closed box welds 36 duplicated corners to 8 nodes', len(wc) // 3 == 8,
      len(wc) // 3)
check('closed box reports no open boundary', M.boundary_loops(wi) == [])

c, i = build(FACES[:2] + FACES[4:])                 # lid removed
wc, wi = M.weld(c, i)
loops = M.boundary_loops(wi)
check('open box finds exactly one boundary loop', len(loops) == 1)

span = M.loop_span_mm(wc, loops[0])
check('gap span is the widest distance across the hole, not the bbox diagonal',
      abs(span - math.sqrt(200)) < 1e-9, span)

fc, fi, filled, left = M.fill_loops(wc, wi, loops, max_gap_mm=20.0, fill_all=False)
check('filling the hole restores a watertight mesh',
      len(filled) == 1 and not M.boundary_loops(fi))

zs = [normal(*[(fc[3 * fi[t + k]], fc[3 * fi[t + k] + 1], fc[3 * fi[t + k] + 2])
               for k in range(3)])[2]
      for t in range(10 * 3, len(fi), 3)]
check('the fill patch winds outward, matching the shell', all(z > 0 for z in zs), zs)

_fc, _fi, filled, left = M.fill_loops(wc, wi, loops, max_gap_mm=5.0, fill_all=False)
check('a hole wider than the limit is deliberately left open',
      not filled and len(left) == 1)

# --------------------------------------------------------------------------- #
# the escalation ladders
# --------------------------------------------------------------------------- #
print('ladders')


def opts(**kw):
    return M.dev_options(kw)


def labels(o, diag=30.0):
    return [r['label'] for r in M.build_ladder(o, diag)]


check('organic ladder climbs accuracy',
      labels(opts(method_key='organic')) == ['Medium', 'High', 'Precise'],
      labels(opts(method_key='organic')))
check('prismatic ladder tightens the face groups',
      len(labels(opts(method_key='prismatic'))) == 3)
check('faceted ladder can drop the reduction',
      labels(opts(method_key='faceted'))[-1] == 'no reduction')
check('form ladder steps the cage DOWN when squareness is off',
      [r['cage'] for r in M.build_ladder(
          opts(method_key='form', cage=24, square_iters=0), 30.0)] == [24, 18, 13],
      [r['cage'] for r in M.build_ladder(
          opts(method_key='form', cage=24, square_iters=0), 30.0)])
# With squareness on the ladder gains a second dimension. Shape is backed off
# BEFORE cage size, because cage size is what carries anatomical fidelity and a
# squareness force is the thing Fusion's T-Spline validator actually refuses.
_sq = M.build_ladder(opts(method_key='form', cage=24, square_iters=120), 30.0)
check('form ladder backs squareness off before it touches the cage',
      [r['cage'] for r in _sq] == [24, 24, 24]
      and [r['square_scale'] for r in _sq] == [1.0, 0.5, 0.25],
      [(r['cage'], r['square_scale']) for r in _sq])
check('no ladder is ever empty',
      all(M.build_ladder(opts(method_key=k), 30.0)
          for k in ('organic', 'prismatic', 'faceted', 'form', 'custom_organic')))
check('a zero-size part does not divide by zero',
      bool(M.build_ladder(opts(method_key='custom_organic'), 0.0)))

o = M.dev_options({'feature_angle_deg': 12.0, 'repair': 'One-touch fix'})
check('degrees convert to radians at the option boundary',
      abs(o['feature_angle_rad'] - math.radians(12)) < 1e-12)
check('repair mode resolves by label', o['repair_type'] is not None)

# --------------------------------------------------------------------------- #
# report rendering, including the degenerate case a one-sided test used to pass
# --------------------------------------------------------------------------- #
print('report')
rep = {
    'source': 'm', 'result': 'm [part]', 'tri_in': 100, 'nodes_in': 52,
    'diag_mm': 1063.0, 'extent_mm': [509.0, 767.7, 530.9], 'closed_in': True,
    'oriented_in': True, 'loops': [], 'filled': [], 'left_open': [],
    'tri_prepared': 100, 'ladder': ['a', 'b'], 'rung_used': 'b',
    'attempts': [
        {'rung': 'a', 'max_mm': None, 'faces': 0, 'error': 'MESH_FAILED_BREP'},
        {'rung': 'b', 'max_mm': 633.4, 'out_mm': 0.0, 'back_mm': 633.4,
         'faces': 1, 'solid': False}],
    'deviation': {
        'out': {'n': 10, 'max': 0.0, 'mean': 0, 'rms': 0, 'p95': 0, 'p99': 0},
        'back': {'n': 10, 'max': 633.4, 'mean': 1, 'rms': 1, 'p95': 1, 'p99': 1}},
    'is_form': False, 'faces_out': 1, 'solid_out': False, 'volume_out_mm3': None,
    'volume_in_mm3': 1.0, 'warnings': [], 'version': 't',
    'methods_available': ['faceted', 'prismatic'],
}
txt = M.format_report(rep, opts(method_key='prismatic', max_error_mm=0.5), 1.0)
check('a body that dropped the whole part is reported OVER budget',
      'OVER the 0.5000 mm budget' in txt)
check('a non-solid result from a closed mesh is called out',
      'dropped geometry' in txt)
check('a refused rung shows its reason', 'refused: MESH_FAILED_BREP' in txt)


# --------------------------------------------------------------------------- #
# cylinder / cone fitting, on tessellated shapes built here (mm in, cm to the fitter)
# --------------------------------------------------------------------------- #

def prism_mesh(n, ra, rb, length, axis, origin, centre_fan=False, phase=0.0):
    ln = math.sqrt(sum(a * a for a in axis))
    ax = [a / ln for a in axis]
    h = [1.0, 0.0, 0.0] if abs(ax[0]) < 0.9 else [0.0, 1.0, 0.0]
    e1 = [ax[1] * h[2] - ax[2] * h[1], ax[2] * h[0] - ax[0] * h[2], ax[0] * h[1] - ax[1] * h[0]]
    l1 = math.sqrt(sum(a * a for a in e1))
    e1 = [a / l1 for a in e1]
    e2 = [ax[1] * e1[2] - ax[2] * e1[1], ax[2] * e1[0] - ax[0] * e1[2], ax[0] * e1[1] - ax[1] * e1[0]]
    V, T = [], []
    for z, r in ((0.0, ra), (length, rb)):
        for k in range(n):
            a = 2 * math.pi * k / n + phase
            V.append([origin[i] + z * ax[i] + r * (math.cos(a) * e1[i] + math.sin(a) * e2[i])
                      for i in range(3)])
    for k in range(n):
        j = (k + 1) % n
        T += [(k, j, n + j), (k, n + j, n + k)]
    if centre_fan:
        V.append(list(origin))
        V.append([origin[i] + length * ax[i] for i in range(3)])
        for k in range(n):
            j = (k + 1) % n
            T += [(2 * n, j, k), (2 * n + 1, n + k, n + j)]
    else:
        for k in range(1, n - 1):
            T += [(0, k + 1, k), (n, n + k, n + k + 1)]
    return [c * 0.1 for v in V for c in v], [i for t in T for i in t]


print('cylinder / cone fit')
c, i = prism_mesh(24, 2.5, 2.5, 30.0, (0.3, 0.5, 0.8), (10.0, -5.0, 3.0))
f = M.fit_cylinder_cone(c, i)
check('a 24-facet cylinder on a skew axis is a cylinder', f['kind'] == 'cylinder')
check('its diameter is 5 mm', abs(2 * f['r_a'] * 10 - 5.0) < 1e-6, 2 * f['r_a'] * 10)
check('its length is 30 mm', abs(f['length'] * 10 - 30.0) < 1e-6, f['length'] * 10)
check('its facets are counted', f['segments'] == 24, f['segments'])
c, i = prism_mesh(21, 1.0, 1.0, 40.0, (0, 0, 1), (0.0, 0.0, 0.0), centre_fan=True)
f = M.fit_cylinder_cone(c, i)
check('caps fanned from a centre vertex still fit', abs(2 * f['r_a'] * 10 - 2.0) < 1e-6)
c, i = prism_mesh(32, 4.0, 2.0, 20.0, (1.0, 0.2, -0.1), (10.0, -5.0, 3.0))
f = M.fit_cylinder_cone(c, i)
ds = sorted((2 * f['r_a'] * 10, 2 * f['r_b'] * 10))   # the axis may point either way
check('a cone is a cone with both diameters',
      f['kind'] == 'cone' and abs(ds[0] - 4.0) < 1e-6 and abs(ds[1] - 8.0) < 1e-6,
      (f['kind'], ds))
c, i = prism_mesh(48, 10.0, 10.0, 2.0, (0, 1, 0), (0.0, 0.0, 0.0))
f = M.fit_cylinder_cone(c, i)
check('a short disc is still a cylinder of the right length',
      abs(f['length'] * 10 - 2.0) < 1e-6 and abs(2 * f['r_a'] * 10 - 20.0) < 1e-6)
c, i = prism_mesh(4, 5.0, 5.0, 10.0, (0, 0, 1), (0.0, 0.0, 0.0), phase=0.3)
try:
    M.fit_cylinder_cone(c, i)
    check('a box is not a cylinder', False)
except ValueError:
    check('a box is not a cylinder', True)
c, i = prism_mesh(24, 2.5, 2.5, 30.0, (0, 0, 1), (0.0, 0.0, 0.0))
c[3 * 5 + 0] += 0.05                       # one vertex 0.5 mm out of place
try:
    M.fit_cylinder_cone(c, i)
    check('a dented cylinder is refused', False)
except ValueError:
    check('a dented cylinder is refused', True)

print()
if FAILED:
    print('%d FAILED: %s' % (len(FAILED), ', '.join(FAILED)))
    raise SystemExit(1)
print('all geometry tests passed')
