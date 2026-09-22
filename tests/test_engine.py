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


# --------------------------------------------------------------------------- #
# a cage that is not all quads
# --------------------------------------------------------------------------- #
print('mixed-arity cage (the gmsh cage carries triangles)')

# a cube with one face split into two triangles: still closed and manifold
CP = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
               [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
CF = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5),
      (3, 7, 6, 2), (0, 4, 7), (0, 7, 3)]

L = M.limit_operator(len(CP), CF)
rows = np.asarray(L.sum(axis=1)).ravel()
check('the mixed-arity limit operator is affine (rows sum to 1)',
      np.allclose(rows, 1.0), '%.3e' % np.abs(rows - 1.0).max())

# repeated subdivision converges to the same limit, and originals keep their index
p_sd, f_sd = CP.copy(), CF
for _ in range(7):
    p_sd, f_sd = M.catmull_clark(p_sd, f_sd)
err = np.abs((L @ CP) - p_sd[:len(CP)]).max()
check('it agrees with 7 rounds of subdivision', err < 1.0e-4, '%.3e' % err)

# the all-quad path must be untouched by the generalisation
QV, QG = M.cube_sphere(6)
QQ = M.cage_quads(QG)
direct = M.limit_matrix(len(QV), M.quad_arrays(QG),
                        M.cage_neighbours(len(QV), QG))
same = abs(direct - M.limit_operator(len(QV), QQ)).max()
check('an all-quad cage still takes the closed-form mask exactly',
      same == 0.0, same)

lim = M.limit_surface_mesh(CP, CF, levels=2)[0]
check('a mixed cage gives a finite limit surface (this is what produced NaN)',
      np.isfinite(lim).all())

# --------------------------------------------------------------------------- #
# the scan-derived quad cage
# --------------------------------------------------------------------------- #
print('scan-derived quad cage')

qv, qf = M.quadrangulate(np.asarray(sphere.vertices),
                         [tuple(int(i) for i in t) for t in sphere.faces])
check('quadrangulating a triangle mesh gives three quads per triangle',
      len(qf) == 3 * len(sphere.faces), len(qf))
check('every face is a quad', all(len(f) == 4 for f in qf))
check('quadrangulation moves nothing: the originals keep their position',
      np.allclose(qv[:len(sphere.vertices)], sphere.vertices))
qtri = np.array([(f[0], f[1], f[2]) for f in qf]
                + [(f[0], f[2], f[3]) for f in qf], np.int64)
check('the quadrangulated cage is watertight',
      trimesh.Trimesh(vertices=qv, faces=qtri, process=False).is_watertight)

# 2 mm on a 30 mm sphere, which is what the dialog now defaults to
res = run(sphere, cage_source='scan', hmax_mm=2.0, output='tsm', limit_iters=15)
check('the scan cage produces a Form output', res['ok'],
      res.get('error', '')[-200:])
if res['ok']:
    s = res['stats']
    check('it reports which cage it used', s.get('cage_source') == 'scan',
          s.get('cage_source'))
    check('it measures how well the cage covers the scan',
          (s.get('cage_cover') or {}).get('p99_mm', 99) < 1.0,
          s.get('cage_cover'))
    check('scan cage max deviation under 0.10 mm',
          s['surf_to_mesh_max_mm'] < 0.10, s['surf_to_mesh_max_mm'])
    check('scan cage has no folded quad', s['cage_folded_quads'] == 0,
          s['cage_folded_quads'])
    check('scan cage volume within 0.1%',
          abs(s['signed_volume_mm3'] - s['mesh_volume_mm3'])
          / s['mesh_volume_mm3'] < 1e-3)
    check('a finer element is asked for and a finer cage arrives',
          s['cage_quads'] > 200, s['cage_quads'])

# the property the gmsh cage never had: asking for a smaller element gives a
# better body. Its density came from a curvature field pinned at the size floor,
# so every rung of the ladder produced the same cage.
coarse = run(sphere, cage_source='scan', hmax_mm=4.0, output='tsm', limit_iters=15)
fine = run(sphere, cage_source='scan', hmax_mm=1.5, output='tsm', limit_iters=15)
if coarse['ok'] and fine['ok']:
    check('a smaller element gives a denser cage',
          fine['stats']['cage_quads'] > coarse['stats']['cage_quads'],
          '%d vs %d' % (fine['stats']['cage_quads'],
                        coarse['stats']['cage_quads']))
    check('and a closer fit',
          fine['stats']['surf_to_mesh_rms_mm']
          < coarse['stats']['surf_to_mesh_rms_mm'],
          '%.4f vs %.4f' % (fine['stats']['surf_to_mesh_rms_mm'],
                            coarse['stats']['surf_to_mesh_rms_mm']))

print()
if FAILED:
    print('%d FAILED: %s' % (len(FAILED), ', '.join(FAILED)))
    raise SystemExit(1)
print('all engine tests passed')
