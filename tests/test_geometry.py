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
check('form ladder steps DOWN, because validity binds before accuracy',
      [r['cage'] for r in M.build_ladder(opts(method_key='form', cage=24), 30.0)]
      == [24, 18, 13],
      [r['cage'] for r in M.build_ladder(opts(method_key='form', cage=24), 30.0)])
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

print()
if FAILED:
    print('%d FAILED: %s' % (len(FAILED), ', '.join(FAILED)))
    raise SystemExit(1)
print('all geometry tests passed')
