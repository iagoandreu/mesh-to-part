"""PROTOTYPE: replace the cube-sphere cage with a gmsh quad mesh.

Not wired into the add-in yet. Run it against a mesh exported with
    py -3.10 tools/mtp.py export_mesh '{"body":"<name>","path":"mesh.json"}'
then feed the TSM to Fusion with tsm_try (keep=true).

Needs gmsh:  pip install gmsh

Measured on a 75k-triangle skin scan against the cube-sphere path:
    element edge ratio   med 1.32  p95 2.44  max  8.95   (was 1.40 / 4.59 / 18.28)
    limit surface        max 0.1938  mean 0.0180 mm      (was max 2.5 to 4.8 mm)
    Fusion body volume   -0.009 %                        (was +0.04 to +0.14 %)
    total time           about 10 s                      (was 70 to 120 s)
"""

import json, math, sys, time
import numpy as np, trimesh, pyvista as pv, gmsh
from collections import defaultdict
sys.path.insert(0, r'c:/Users/iandreu/Documents/Claude_iago/FusionAddIns/MeshToPart')
import mtp_organic as mo

HMAX = float(sys.argv[1]) if len(sys.argv) > 1 else 8.0
HMIN = float(sys.argv[2]) if len(sys.argv) > 2 else 1.5
CURV = float(sys.argv[3]) if len(sys.argv) > 3 else 16.0
OUT  = sys.argv[4] if len(sys.argv) > 4 else 'piel_fitted.tsm'
def log(m): print(m, flush=True)

d = json.load(open('piel.json'))
full = trimesh.Trimesh(vertices=np.asarray(d['coords'],float).reshape(-1,3)*10.0,
                       faces=np.asarray(d['idx'],np.int64).reshape(-1,3), process=True)
fc = np.asarray(full.faces)
pdv = pv.PolyData(np.asarray(full.vertices,float), np.hstack([np.full((len(fc),1),3),fc]).ravel())
sm = pdv.decimate(1.0-20000.0/len(fc)).clean().triangulate()
sfa = np.asarray(sm.faces).reshape(-1,4)[:,1:]
trimesh.Trimesh(vertices=np.asarray(sm.points,float), faces=sfa, process=True).export('piel_full.stl')

gmsh.initialize(); gmsh.option.setNumber('General.Terminal',0)
gmsh.merge('piel_full.stl')
gmsh.model.mesh.classifySurfaces(math.pi, True, True, math.pi)
gmsh.model.mesh.createGeometry(); gmsh.model.geo.synchronize()
for k,v in (('Mesh.MeshSizeMin',HMIN),('Mesh.MeshSizeMax',HMAX),
            ('Mesh.MeshSizeFromCurvature',CURV),('Mesh.MeshSizeExtendFromBoundary',0),
            ('Mesh.Algorithm',8),('Mesh.RecombineAll',1),
            ('Mesh.RecombinationAlgorithm',3),('Mesh.SubdivisionAlgorithm',0)):
    gmsh.option.setNumber(k,v)
gmsh.model.mesh.generate(2)
nt,nc,_ = gmsh.model.mesh.getNodes(); pts = np.asarray(nc,float).reshape(-1,3)
remap = {int(t):i for i,t in enumerate(nt)}
et,_ei,en = gmsh.model.mesh.getElements(2); faces=[]
for t,nodes in zip(et,en):
    k = {2:3,3:4}.get(int(t))
    if k: faces.extend(tuple(int(x) for x in r) for r in
                       np.array([remap[int(x)] for x in nodes],np.int64).reshape(-1,k))
gmsh.finalize()
used = sorted({v for f in faces for v in f}); ix = {v:i for i,v in enumerate(used)}
pts = pts[used]; faces = [tuple(ix[v] for v in f) for f in faces]

exec(open('gq2.py').read().split('def orient(faces):')[1].split('faces = orient(faces)')[0]
     .join(['def orient(faces):','']))
faces = orient(faces)
tri = []
for f in faces:
    tri.append((f[0],f[1],f[2]))
    if len(f)==4: tri.append((f[0],f[2],f[3]))
if trimesh.Trimesh(vertices=pts, faces=np.array(tri), process=False).volume < 0:
    faces = [tuple(reversed(f)) for f in faces]
nq = sum(1 for f in faces if len(f)==4)
log('cage: %d quads + %d tris, %d verts' % (nq, len(faces)-nq, len(pts)))

# --- fit the LIMIT surface to the mesh -------------------------------------
pts_cm = pts * 0.1                                    # engine works in cm
proj = mo.Projector(trimesh.Trimesh(vertices=full.vertices*0.1,
                                    faces=full.faces, process=True))
table = mo.neighbours_from_quads(len(pts_cm), faces)
qs = np.array([f for f in faces if len(f)==4], np.int64)
qa = (qs[:,0].copy(), qs[:,1].copy(), qs[:,2].copy(), qs[:,3].copy())
log('table width %d (max valence), %d quad faces drive the limit operator'
    % (table.shape[1], len(qs)))
t0=time.time()
fit, mx, mn, folds, lam = mo.fit_limit_unfolded(pts_cm, qa, table, proj, 20, 1.0e-3, 60, log)
log('limit fit: surface sits max %.4f mean %.4f mm from the mesh, folds %d, lam %.1e, %.0fs'
    % (mx*10.0, mn*10.0, folds, lam, time.time()-t0))

lim_pts, lim_tris = mo.limit_surface_mesh(fit, [list(f) for f in faces if len(f)==4], levels=2)
ls = trimesh.Trimesh(vertices=np.asarray(lim_pts)*10.0, faces=np.asarray(lim_tris), process=True)
log('limit surface volume %.0f mm3 vs mesh %.0f (%+.2f%%)'
    % (abs(ls.volume), full.volume, 100.0*(abs(ls.volume)-full.volume)/full.volume))

tsm = mo.write_tsm(fit, faces)
open(OUT,'w',encoding='utf-8').write(tsm)
log('TSM %s: %d bytes' % (OUT, len(tsm)))
