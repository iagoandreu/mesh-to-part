"""Engine tests: the fitting maths, on shapes with a known answer.

    py -3.10 tests/test_engine.py

Needs numpy, scipy and trimesh. Does not need Fusion. The point of these is that
the reference figures are analytic, not "whatever it did last time".
"""

import json
import os
import subprocess
import sys
import tempfile

import numpy as np
import trimesh

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

import mtp_organic as M                                     # noqa: E402

ENGINE = os.path.join(ROOT, 'mtp_organic.py')
FAILED = []


def check(name, cond, detail=''):
    print('  %s %s%s' % ('ok  ' if cond else 'FAIL', name,
                         '' if cond else '   %s' % (detail,)))
    if not cond:
        FAILED.append(name)


def run(mesh, **kw):
    job = {'coords': mesh.vertices.ravel().tolist(),
           'idx': mesh.faces.ravel().tolist()}
    job.update(kw)
    d = tempfile.mkdtemp()
    jp, op = os.path.join(d, 'j.json'), os.path.join(d, 'o.json')
    with open(jp, 'w') as fh:
        json.dump(job, fh)
    subprocess.run([sys.executable, ENGINE, jp, op], capture_output=True)
    with open(op) as fh:
        return json.load(fh)


# --------------------------------------------------------------------------- #
# the cage
# --------------------------------------------------------------------------- #
print('cage')
verts, grids = M.cube_sphere(2)
quads = M.cage_quads(grids)
check('a cube-sphere at n=2 has quad-ball topology (26/48/24)',
      (len(verts), len(quads)) == (26, 24), (len(verts), len(quads)))
check('every cage vertex is on the unit sphere',
      abs(np.linalg.norm(verts, axis=1) - 1.0).max() < 1e-12)

# --------------------------------------------------------------------------- #
# TSM: structure must match what Fusion's own Quad Ball uses
# --------------------------------------------------------------------------- #
print('TSM writer')
tsm = M.write_tsm(verts * 2.0, quads)
lines = tsm.splitlines()
L = np.array([[int(x) for x in ln.split()[1:]] for ln in lines if ln.startswith('l ')])
a, b, c, vert, face, edge, _flag = [L[:, i] for i in range(7)]
n = len(L)
nf = len([1 for ln in lines if ln.startswith('f ')])
ne = len([1 for ln in lines if ln.startswith('e ')])
nv = len([1 for ln in lines if ln.startswith('v ')])
check('header is #TS0200', lines[0] == '#TS0200')
check('counts are 24 faces / 48 edges / 26 vertices / 96 links',
      (nf, ne, nv, n) == (24, 48, 26, 96), (nf, ne, nv, n))
check('Euler characteristic is 2', nv - ne + nf == 2)
check('next and prev are inverses', bool((b[a] == np.arange(n)).all()))
check('mate is an involution', bool((c[c] == np.arange(n)).all()))
check('mate shares the edge', bool((edge[c] == edge).all()))
check('mate belongs to the other face', bool((face[c] != face).all()))
check('next stays within the face', bool((face[a] == face).all()))
check('vertex field is the link DESTINATION', bool((vert[c] == vert[b]).all()))

# --------------------------------------------------------------------------- #
# Catmull-Clark limit surface: checked against a figure Fusion itself produced
# --------------------------------------------------------------------------- #
print('limit surface')
table = M.cage_neighbours(len(verts), grids)
qa = M.quad_arrays(grids)
lim = M.limit_positions(verts * 2.0, qa, table)
r = np.linalg.norm(lim, axis=1)
check('limit points of a radius-2 cage land at 1.6457 (Fusion built a 32.9146 mm '
      'body from this cage)', abs(r.min() - 1.6457) < 1e-3, r.min())

lim_pts, lim_tris = M.limit_surface_mesh(verts * 2.0, quads, levels=2)
size = 2 * np.abs(lim_pts).max()
check('subdivided limit surface predicts 32.9146 mm', abs(size * 10 - 32.9146) < 1e-3,
      size * 10)

# --------------------------------------------------------------------------- #
# end to end, both outputs, against analytic shapes
# --------------------------------------------------------------------------- #
print('fitting a sphere (radius 15 mm)')
sphere = trimesh.creation.icosphere(subdivisions=4, radius=1.5)

res = run(sphere, cage=16, wrap_iters=50, relax_iters=50, ctrl=17, fair=1e-2)
check('NURBS output succeeds', res['ok'], res.get('error', '')[-200:])
if res['ok']:
    s = res['stats']
    check('NURBS max deviation under 0.05 mm', s['surf_to_mesh_max_mm'] < 0.05,
          s['surf_to_mesh_max_mm'])
    check('NURBS rms under 0.01 mm', s['surf_to_mesh_rms_mm'] < 0.01,
          s['surf_to_mesh_rms_mm'])
    check('NURBS cage has no folded quad', s['cage_folded_quads'] == 0,
          s['cage_folded_quads'])
    check('NURBS volume within 0.1%',
          abs(s['signed_volume_mm3'] - s['mesh_volume_mm3'])
          / s['mesh_volume_mm3'] < 1e-3)

res = run(sphere, cage=16, wrap_iters=50, relax_iters=50, output='tsm',
          limit_iters=15)
check('Form output succeeds', res['ok'], res.get('error', '')[-200:])
if res['ok']:
    s = res['stats']
    check('Form max deviation under 0.10 mm', s['surf_to_mesh_max_mm'] < 0.10,
          s['surf_to_mesh_max_mm'])
    check('Form cage has no folded quad, or Fusion would refuse it',
          s['cage_folded_quads'] == 0, s['cage_folded_quads'])
    check('Form volume within 0.1%',
          abs(s['signed_volume_mm3'] - s['mesh_volume_mm3'])
          / s['mesh_volume_mm3'] < 1e-3)
    check('the emitted TSM starts with the right magic',
          res['tsm'].startswith('#TS0200'))

print()
if FAILED:
    print('%d FAILED: %s' % (len(FAILED), ', '.join(FAILED)))
    raise SystemExit(1)
print('all engine tests passed')
