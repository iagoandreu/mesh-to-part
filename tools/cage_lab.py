"""Compare T-spline control cages on a scan, each fitted by the engine and measured.

    py -3.10 tools/cage_lab.py scan.stl <element_mm> scan instant pair dominant

scan / instant are the engine's cages; pair / dominant are MeshLab quad remeshes
(pip install pymeshlab), kept for comparison. Reports the fit deviation both ways,
the extraordinary-vertex share, folds and the limit surface's normal noise.
"""
import os
import sys
import time
import json

import numpy as np
import trimesh

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import mtp_organic as mo


def meshlab_cage(mesh_cm, h_mm, method):
    """Isotropic remesh at h, then MeshLab's triangle pairing or quad-dominant pass."""
    import pymeshlab
    ms = pymeshlab.MeshSet()
    ms.add_mesh(pymeshlab.Mesh(np.asarray(mesh_cm.vertices) * 10.0, np.asarray(mesh_cm.faces)))
    ms.meshing_isotropic_explicit_remeshing(targetlen=pymeshlab.PureValue(h_mm),
                                            iterations=12, featuredeg=60.0)
    if method == 'pair':
        ms.meshing_tri_to_quad_by_smart_triangle_pairing()
    else:
        ms.meshing_tri_to_quad_dominant(level=1)
    m = ms.current_mesh()
    pts = m.vertex_matrix()
    faces = mo.orient_faces([tuple(int(i) for i in p) for p in m.polygonal_face_list()])
    tri = [(q[0], q[k - 1], q[k]) for q in faces for k in range(2, len(q))]
    if trimesh.Trimesh(vertices=pts, faces=np.array(tri), process=False).volume < 0:
        faces = [tuple(reversed(q)) for q in faces]
    return pts * 0.1, faces


def normal_noise(v, f):
    """p95/p99 angle (deg) between each vertex normal and its 1-ring's mean normal."""
    tm = trimesh.Trimesh(vertices=v, faces=f, process=False)
    n = tm.vertex_normals
    e = tm.edges_unique
    acc = np.zeros_like(n)
    np.add.at(acc, e[:, 0], n[e[:, 1]])
    np.add.at(acc, e[:, 1], n[e[:, 0]])
    m = acc / np.linalg.norm(acc, axis=1, keepdims=True)
    a = np.degrees(np.arccos(np.clip((n * m).sum(1), -1, 1)))
    return round(float(np.percentile(a, 95)), 2), round(float(np.percentile(a, 99)), 2)


def main():
    path, h = sys.argv[1], float(sys.argv[2])
    raw = trimesh.load(path)
    mesh = trimesh.Trimesh(vertices=raw.vertices * 0.1, faces=raw.faces, process=True)
    print('scan %d tris, watertight %s, area %.0f mm2'
          % (len(mesh.faces), mesh.is_watertight, mesh.area * 100))
    for meth in sys.argv[3:]:
        t0 = time.time()
        job = {'hmax_mm': h, 'output': 'tsm', 'cage_source': 'scan',
               'limit_iters': 20, 'limit_fair': 1e-3,
               'coords': mesh.vertices.ravel().tolist(), 'idx': mesh.faces.ravel().tolist()}
        if meth == 'scan':
            pts, faces = mo.scan_quad_cage(mesh, job, lambda s: None)
        elif meth == 'instant':
            pts, faces = mo.instant_quad_cage(mesh, job, lambda s: None)
        else:
            pts, faces = meshlab_cage(mesh, h, meth)
        val = np.bincount(np.concatenate([np.array(q) for q in faces]), minlength=len(pts))
        rep = {'verts': len(pts), 'faces': len(faces),
               'extraordinary_pct': round(float((val != 4).mean() * 100), 1)}
        orig = mo.scan_quad_cage
        mo.scan_quad_cage = lambda m, j, l, _p=pts, _f=faces: (_p, _f)
        try:
            res = mo.build(dict(job), lambda s: None)
        finally:
            mo.scan_quad_cage = orig
        st = res['stats']
        for k in ('cage_folded_quads', 'surf_to_mesh_max_mm', 'mesh_to_surf_max_mm',
                  'mesh_to_surf_p99_mm', 'mesh_to_surf_rms_mm'):
            rep[k] = round(st[k], 3) if isinstance(st[k], float) else st[k]
        lp, lt = mo.limit_surface_mesh(np.array(res['cage_points']), faces, 2)
        rep['noise_p95_p99_deg'] = normal_noise(lp * 10, lt)
        rep['secs'] = round(time.time() - t0, 1)
        print(meth, json.dumps(rep), flush=True)


if __name__ == '__main__':
    main()
