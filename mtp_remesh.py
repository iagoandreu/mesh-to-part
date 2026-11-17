"""Mesh to Part - Organic Form by quad remesh. The fitting engine.

Runs in the external Python (numpy, scipy, trimesh, open3d, PyNanoInstantMeshes):

    python mtp_remesh.py job.json out.json

The one organic method. A field-aligned quad cage from Instant Meshes (about 3 %
extraordinary vertices, so the Form converts to a few hundred B-Rep faces, not
one per element) and its Catmull-Clark limit surface fitted to the mesh by dense
point-to-plane least squares in both directions. Everything after loading costs
in proportion to the cage, not to the input triangles: a 6.5 k and a 173 k
triangle femur both take about 10 s.

job['mode'] == 'measure' compares a tessellation of the finished body against the
source instead, so the report states what Fusion actually built.

Progress lines go to job['progress'] (one line, overwritten) for the add-in's
progress dialog; the add-in kills this process when the user cancels.
"""

import json
import os
import sys
import time

import numpy as np
import trimesh
from scipy.sparse import csr_matrix, diags, hstack, vstack
from scipy.sparse.linalg import LinearOperator, cg
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mtp_organic as O          # subdivision / limit operators, TSM writer

PROGRESS = None


def progress(msg):
    if PROGRESS:
        try:
            with open(PROGRESS + '.part', 'w', encoding='utf-8') as fh:
                fh.write(msg)
            os.replace(PROGRESS + '.part', PROGRESS)
        except OSError:
            pass


def quiet(fn, *a, **k):
    """Instant Meshes prints its whole log on the C stdout; keep it off ours."""
    sys.stdout.flush()
    fd = os.dup(1)
    nul = os.open(os.devnull, os.O_WRONLY)
    os.dup2(nul, 1)
    try:
        return fn(*a, **k)
    finally:
        os.dup2(fd, 1)
        os.close(nul)
        os.close(fd)


def load(vpath, fpath):
    v = np.fromfile(vpath, np.float64).reshape(-1, 3)
    f = np.fromfile(fpath, np.int32).reshape(-1, 3)
    return trimesh.Trimesh(v, f, process=True)


class Surface(object):
    """Exact closest points on a triangle mesh (open3d's BVH)."""

    def __init__(self, mesh):
        import open3d as o3d
        self.o3d = o3d
        self.scene = o3d.t.geometry.RaycastingScene()
        self.scene.add_triangles(o3d.core.Tensor(np.asarray(mesh.vertices, np.float32)),
                                 o3d.core.Tensor(np.asarray(mesh.faces, np.uint32)))
        self.normals = np.asarray(mesh.face_normals, float)

    def closest(self, X):
        r = self.scene.compute_closest_points(self.o3d.core.Tensor(np.asarray(X, np.float32)))
        T = r['points'].numpy().astype(float)
        fi = r['primitive_ids'].numpy().astype(np.int64)
        return T, np.linalg.norm(T - X, axis=1), fi

    def distance(self, X):
        return self.scene.compute_distance(
            self.o3d.core.Tensor(np.asarray(X, np.float32))).numpy().astype(float)


# --------------------------------------------------------------------------- #
# cage
# --------------------------------------------------------------------------- #

def quad_cage(mesh, quads):
    """Instant Meshes, about `quads` quads. Its extraction ends in one regular
    subdivision, so the vertex count asked for is a quarter of that."""
    import pynanoinstantmeshes as pim
    V, F = quiet(pim.remesh, np.asarray(mesh.vertices, np.float32),
                 np.asarray(mesh.faces, np.uint32), max(50, int(quads) // 4),
                 rosy=4, posy=4, smooth_iter=2, deterministic=True)
    V = np.asarray(V, float)
    faces = []
    for q in np.asarray(F):
        u = []
        for i in (int(x) for x in q):
            if not u or u[-1] != i:
                u.append(i)
        if len(u) > 1 and u[0] == u[-1]:
            u.pop()
        if len(u) >= 3:
            faces.append(tuple(u))
    faces, _ = O.dissolve_doublets(faces)
    used = np.unique(np.concatenate([np.array(q) for q in faces]))
    remap = -np.ones(len(V), np.int64)
    remap[used] = np.arange(len(used))
    V = V[used]
    faces = [tuple(int(remap[i]) for i in q) for q in faces]
    uses = {}
    for q in faces:
        for k in range(len(q)):
            e = (min(q[k], q[k - 1]), max(q[k], q[k - 1]))
            uses[e] = uses.get(e, 0) + 1
    if any(c != 2 for c in uses.values()):
        raise RuntimeError('the quad remesh came back open; is the mesh closed?')
    faces = O.orient_faces(faces)
    tri = [(q[0], q[k - 1], q[k]) for q in faces for k in range(2, len(q))]
    if trimesh.Trimesh(V, np.array(tri), process=False).volume < 0:
        faces = [tuple(reversed(q)) for q in faces]
    return V, faces


def eval_operator(n, faces, levels):
    """Rows over the cage giving the limit points of the cage subdivided `levels`
    times - a dense, exact sample of the surface Fusion builds - and the
    triangles that connect them (outward)."""
    S = None
    f = [tuple(q) for q in faces]
    nn = n
    for _ in range(levels):
        Si, f = O.subdivision_matrix(nn, f)
        nn = Si.shape[0]
        S = Si if S is None else Si @ S
    tab = O.neighbours_from_quads(nn, f)
    qa = tuple(np.array([x[k] for x in f], np.int64) for k in range(4))
    L = O.limit_matrix(nn, qa, tab)
    E = (L @ S).tocsr() if S is not None else L.tocsr()
    fa = np.array(f, np.int64)
    tris = np.vstack([fa[:, [0, 1, 2]], fa[:, [0, 2, 3]]])
    return E, tris


def block3(M):
    Z = csr_matrix(M.shape)
    return vstack([hstack([M, Z, Z]), hstack([Z, M, Z]), hstack([Z, Z, M])]).tocsr()


def solve(M, r, x0):
    dg = M.diagonal()
    P = LinearOperator(M.shape, matvec=lambda x: x / dg)
    x, _info = cg(M, r, x0=x0, rtol=1e-7, maxiter=500, M=P)
    return x


def fit(mesh, surf, V, faces, iters=6, mu=1e-3, alpha=0.05, beta=1.0, samples=60000,
        log=None, tag=''):
    """Control cage whose limit surface lies on the mesh.

    Point-to-plane both ways: every surface sample to its closest mesh point,
    and every mesh sample to its nearest surface sample - the second term is
    what stops the surface cutting across a ridge. A weak point-to-point term
    holds the samples from sliding, and a graph Laplacian on the control net
    keeps it tame."""
    n = len(V)
    E, _tris = eval_operator(n, faces, 2)
    D = O.graph_laplacian(n, O.neighbours_from_quads(n, faces))
    Db = block3(D)
    Y, fy = trimesh.sample.sample_surface(mesh, samples, seed=1)
    Ny = surf.normals[fy]
    Eb = block3(E)
    C = np.asarray(V, float).copy()
    kd_iters = []
    for it in range(iters):
        progress('Fitting the surface%s (%d/%d)' % (tag, it + 1, iters))
        X = E @ C
        T, d, fi = surf.closest(X)
        Nt = surf.normals[fi]
        A1 = hstack([diags(Nt[:, 0]) @ E, diags(Nt[:, 1]) @ E, diags(Nt[:, 2]) @ E]).tocsr()
        b1 = (Nt * T).sum(1)
        _dd, j = cKDTree(X).query(Y, workers=-1)
        Ej = E[j]
        A2 = hstack([diags(Ny[:, 0]) @ Ej, diags(Ny[:, 1]) @ Ej, diags(Ny[:, 2]) @ Ej]).tocsr()
        b2 = (Ny * Y).sum(1)
        bT = np.concatenate([T[:, 0], T[:, 1], T[:, 2]])
        w1, w2, w3 = 1.0 / len(X), beta / len(Y), alpha / len(X)
        M = (w1 * (A1.T @ A1) + w2 * (A2.T @ A2) + w3 * (Eb.T @ Eb)
             + (mu / n) * (Db.T @ Db)).tocsr()
        r = w1 * (A1.T @ b1) + w2 * (A2.T @ b2) + w3 * (Eb.T @ bT)
        c = solve(M, r, np.concatenate([C[:, 0], C[:, 1], C[:, 2]]))
        C = np.column_stack([c[:n], c[n:2 * n], c[2 * n:]])
        kd_iters.append(float(d.mean()))
        if log:
            log('  fit %d/%d: surface to mesh mean %.4f max %.3f mm'
                % (it + 1, iters, d.mean(), d.max()))
    return C


def offset_cage(C, faces, d, mu=1e-3):
    """The same cage refitted so its limit surface is the exact normal offset of
    the fitted one by d mm: correspondence is by parameter, so it is one linear
    solve per coordinate."""
    n = len(C)
    E, tris = eval_operator(n, faces, 2)
    X = E @ C
    N = trimesh.Trimesh(X, tris, process=False).vertex_normals
    T = X + d * N
    D = O.graph_laplacian(n, O.neighbours_from_quads(n, faces))
    M = (E.T @ E / len(X) + (mu / n) * (D.T @ D)).tocsr()
    out = np.empty_like(C)
    for k in range(3):
        out[:, k] = solve(M, E.T @ T[:, k] / len(X), C[:, k])
    return out


# --------------------------------------------------------------------------- #
# measuring
# --------------------------------------------------------------------------- #

def stats(d):
    return {'max': float(d.max()), 'p99': float(np.percentile(d, 99)),
            'p95': float(np.percentile(d, 95)), 'rms': float(np.sqrt((d ** 2).mean()))}


def worst_spots(P, d, above, count=4):
    """Where the largest deviations are: the worst point of each cluster of points
    over `above`, clusters 5 mm apart, worst first."""
    k = d > above
    if not k.any():
        return []
    Q, dq = P[k], d[k]
    order = np.argsort(-dq)
    spots = []
    for i in order:
        if all(np.linalg.norm(Q[i] - s[0]) > 5.0 for s in spots):
            spots.append((Q[i], dq[i], int(((np.linalg.norm(Q - Q[i], axis=1) < 5.0)).sum())))
            if len(spots) >= count:
                break
    return [{'at_mm': [round(float(x), 2) for x in p], 'dev_mm': round(float(v), 3),
             'points': c} for p, v, c in spots]


def compare(part, mesh, samples=200000):
    """Both directions between two triangle meshes: part -> mesh and mesh -> part,
    over area-uniform samples plus the vertices."""
    pa = np.vstack([trimesh.sample.sample_surface(part, samples, seed=3)[0], part.vertices])
    pb = np.vstack([trimesh.sample.sample_surface(mesh, samples, seed=4)[0], mesh.vertices])
    d1 = Surface(mesh).distance(pa)
    d2 = Surface(part).distance(pb)
    s1, s2 = stats(d1), stats(d2)
    worst = max(s1['p99'], s2['p99'])
    return {'part_to_mesh': s1, 'mesh_to_part': s2,
            'worst_p99': worst, 'worst_max': max(s1['max'], s2['max']),
            'spots': worst_spots(np.vstack([pa, pb]), np.concatenate([d1, d2]),
                                 max(0.3, 3.0 * worst))}


def limit_mesh(C, faces, levels=3):
    E, tris = eval_operator(len(C), faces, levels)
    return trimesh.Trimesh(E @ C, tris, process=False)


# --------------------------------------------------------------------------- #
# the job
# --------------------------------------------------------------------------- #

def auto_element(area_mm2, quads=8000):
    """Element edge for about `quads` quads, kept between 0.8 and 6 mm."""
    return float(np.clip(np.sqrt(area_mm2 / quads), 0.8, 6.0))


def run(job, log):
    t0 = time.time()
    progress('Reading the mesh')
    mesh = load(job['v'], job['f'])
    if not mesh.is_watertight and job.get('fill_gaps', True):
        trimesh.repair.fill_holes(mesh)
        log('open boundaries filled')
    if not mesh.is_watertight:
        raise RuntimeError('the mesh is not closed (watertight) and its holes could not '
                           'be filled; repair it first (Mesh > Repair)')
    surf = Surface(mesh)
    area = float(mesh.area)
    h = float(job.get('element_mm') or 0.0) or auto_element(area)
    budget = float(job.get('budget_mm') or 0.0)
    limit_s = float(job.get('time_s', 60.0))
    log('mesh %d triangles, %.0f mm2, element %.2f mm' % (len(mesh.faces), area, h))
    tries = []
    best = None
    while True:
        t1 = time.time()
        progress('Quad remesh, element %.2f mm' % h)
        V, faces = quad_cage(mesh, area / (h * h))
        val = np.bincount(np.concatenate([np.array(q) for q in faces]), minlength=len(V))
        star = int((val != 4).sum())
        log('cage: %d quads, %d vertices, %d extraordinary (%.1f %%)'
            % (len(faces), len(V), star, 100.0 * star / len(V)))
        C = fit(mesh, surf, V, faces, log=log, tag=', element %.2f mm' % h)
        progress('Measuring, element %.2f mm' % h)
        cmp_ = compare(limit_mesh(C, faces), mesh)
        dt = time.time() - t1
        row = {'element_mm': round(h, 3), 'quads': len(faces), 'extraordinary': star,
               'p99_mm': round(cmp_['worst_p99'], 4), 'max_mm': round(cmp_['worst_max'], 3),
               'secs': round(dt, 1)}
        tries.append(row)
        log('element %.2f mm: p99 %.4f max %.3f mm in %.1f s'
            % (h, cmp_['worst_p99'], cmp_['worst_max'], dt))
        if best is None or cmp_['worst_p99'] < best[2]['worst_p99']:
            best = (C, faces, cmp_, h, star)
        spent = time.time() - t0
        # a finer cage costs about twice this one; only if it fits in the time
        if (not budget or cmp_['worst_p99'] <= budget or h <= 0.8
                or spent + 2.2 * dt > limit_s):
            break
        h = max(0.8, h * 0.72)

    C, faces, cmp_, h, star = best
    offset = float(job.get('offset_mm') or 0.0)
    off_rep = None
    if offset:
        progress('Offsetting the surface by %+.3f mm' % offset)
        base = limit_mesh(C, faces)
        C = offset_cage(C, faces, offset)
        moved = limit_mesh(C, faces)
        dd = Surface(base).distance(moved.vertices)
        off_rep = {'asked_mm': offset, 'mean_mm': float(dd.mean()), 'min_mm': float(dd.min()),
                   'max_mm': float(dd.max())}
        log('offset %+.3f mm: achieved mean %.4f min %.4f max %.4f mm'
            % (offset, dd.mean(), dd.min(), dd.max()))
    qa = tuple(np.array([q[k] for q in faces], np.int64) for k in range(4))
    folded, _l = O.folded_quads(C, qa)
    with open(job['tsm'], 'w', encoding='utf-8') as fh:
        fh.write(O.write_tsm(C * 0.1, faces))
    lm = limit_mesh(C, faces, 2)
    return {'tsm': job['tsm'], 'element_mm': h, 'quads': len(faces), 'cage_vertices': len(C),
            'extraordinary': star, 'folded_quads': int(folded.sum()), 'tries': tries,
            'engine_check': cmp_, 'offset': off_rep, 'triangles_in': int(len(mesh.faces)),
            'mesh_volume_mm3': float(mesh.volume), 'surface_volume_mm3': float(lm.volume),
            'secs': round(time.time() - t0, 1)}


def run_measure(job, log):
    progress('Measuring the body against the mesh')
    part = load(job['part_v'], job['part_f'])
    mesh = load(job['v'], job['f'])
    out = compare(part, mesh)
    out['part_volume_mm3'] = float(part.volume) if part.is_watertight else None
    out['mesh_volume_mm3'] = float(mesh.volume)
    return out


def main():
    global PROGRESS
    job = json.load(open(sys.argv[1], encoding='utf-8'))
    PROGRESS = job.get('progress')
    lines = []

    def log(msg):
        lines.append(msg)

    try:
        res = run_measure(job, log) if job.get('mode') == 'measure' else run(job, log)
        res['ok'] = True
    except Exception:
        import traceback
        res = {'ok': False, 'error': traceback.format_exc()}
    res['log'] = lines
    with open(sys.argv[2] + '.part', 'w', encoding='utf-8') as fh:
        json.dump(res, fh)
    os.replace(sys.argv[2] + '.part', sys.argv[2])


if __name__ == '__main__':
    main()
