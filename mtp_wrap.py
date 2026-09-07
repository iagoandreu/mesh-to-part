"""
Offset ("wrap") engine for Mesh to Part. Runs OUTSIDE Fusion, under a Python
with numpy, scipy and pyvista -- Fusion's own interpreter has none of them.

WHY THIS EXISTS

Fusion's Offset Face and Shell work in BRep space: they offset each face and
then try to heal the result. On an organic body that fails, because an outward
offset of a concave region self-intersects -- between toes, behind a malleolus,
anywhere the local concave radius is smaller than the offset distance. There is
no surface for Fusion to keep there; the correct answer is that the material
merges, and a per-face offset cannot express that.

A level set can. The signed distance field of a closed surface is defined
everywhere, and the isosurface at distance d is the exact outward offset --
self-merging in concavities is automatic, not a special case. So:

    surface -> voxel occupancy -> signed distance -> isosurface at d -> smooth

The result cannot self-intersect by construction, which is the whole point.

ACCURACY

The distance transform measures to voxel CENTRES, so the field carries a
quantisation of about pitch/2 near the surface. Taubin smoothing removes the
resulting stair pattern without the shrinkage a Laplacian would cause. The
engine measures what it actually achieved -- the real distance from the offset
surface back to the source -- and reports it, rather than assuming it hit d.

Job (JSON):
    coords        flat [x,y,z,...] in CENTIMETRES, as Fusion hands them over
    idx           flat triangle indices
    offset_mm     outward offset distance (negative shrinks)
    pitch_mm      voxel size; 0 or absent picks one from the voxel budget
    close_mm      morphological gap closing before the offset (0 = off)
    smooth_iters  Taubin passes on the extracted isosurface
    budget        maximum voxel count, default 15e6

Result (JSON): coords (CENTIMETRES) + idx for the offset mesh, plus stats.
"""

import json
import sys
import time

import numpy as np

CM = 10.0          # cm -> mm
MM = 0.1           # mm -> cm
DEFAULT_BUDGET = 15.0e6

_LOG = []


def log(msg):
    _LOG.append(msg)
    # ASCII only: this console is cp1252 and a stray arrow kills the run
    print(msg)
    sys.stdout.flush()


def pick_pitch(size_mm, offset_mm, pitch_mm, budget):
    """Voxel size: fine enough to resolve the offset, coarse enough to fit.

    The offset has to span several voxels or the isosurface is a staircase, so
    the ideal is offset/4. That is capped by the voxel budget, because the grid
    is the whole bounding box plus a margin and memory goes as pitch^-3.
    """
    span = size_mm + 2.0 * (abs(offset_mm) + 4.0 * max(pitch_mm, 0.5))
    ideal = max(abs(offset_mm) / 4.0, 1.0e-6) if offset_mm else np.max(size_mm) / 300.0
    if pitch_mm and pitch_mm > 0:
        p = float(pitch_mm)
    else:
        p = float(np.clip(ideal, np.max(size_mm) / 600.0, np.max(size_mm) / 60.0))
    while np.prod(np.ceil(span / p) + 1) > budget:
        p *= 1.15
    return p


def occupancy(mesh, origin, dims, pitch):
    """Solid voxel mask of a watertight mesh.

    trimesh's own voxelizer subdivides triangles until they are smaller than a
    voxel, which is exact for the shell but leaves the interior empty; the fill
    is a flood from outside the padded grid, so an interior cavity stays a
    cavity instead of being swallowed.
    """
    import trimesh
    from scipy import ndimage

    shell = np.zeros(dims, dtype=bool)
    v = mesh.vertices
    f = mesh.faces
    # subdivide long triangles so no voxel is skipped between vertices
    verts, faces = trimesh.remesh.subdivide_to_size(v, f, max_edge=pitch * 0.5)
    ijk = np.floor((verts - origin) / pitch).astype(np.int64)
    np.clip(ijk, 0, np.array(dims) - 1, out=ijk)
    shell[ijk[:, 0], ijk[:, 1], ijk[:, 2]] = True

    outside = ~shell
    labels, _n = ndimage.label(outside)
    # the grid is padded, so index [0,0,0] is guaranteed to be outside the part
    outer = labels[0, 0, 0]
    solid = ~((labels == outer) & outside)
    return solid


def signed_distance(solid, pitch):
    """Positive outside, negative inside, in millimetres."""
    from scipy import ndimage
    d_out = ndimage.distance_transform_edt(~solid, sampling=pitch)
    d_in = ndimage.distance_transform_edt(solid, sampling=pitch)
    return (d_out - d_in).astype(np.float32)


def close_gaps(solid, pitch, radius_mm):
    """Morphological closing expressed on the level set.

    dilate(X, r) is {sdf <= r} and erode(Y, r) is {sdf_Y <= -r}, so a closing is
    two distance transforms rather than a structuring element -- which matters,
    because a ball element of radius r/pitch voxels gets expensive fast.
    """
    grown = signed_distance(solid, pitch) <= radius_mm
    return signed_distance(grown, pitch) <= -radius_mm


def as_polydata(mesh):
    import pyvista as pv
    f = np.asarray(mesh.faces)
    return pv.PolyData(np.asarray(mesh.vertices, float),
                       np.hstack([np.full((len(f), 1), 3), f]).ravel())


def exact_signed_distance(points, mesh):
    """True signed distance from each point to the surface, negative inside.

    trimesh's own closest_point needs rtree, which is not part of this add-in's
    stated dependency set (numpy, scipy, trimesh). VTK ships the same thing in
    vtkImplicitPolyDataDistance and pyvista is already required for the
    marching cubes, so this costs no new dependency.
    """
    import pyvista as pv
    cloud = pv.PolyData(np.asarray(points, float))
    out = cloud.compute_implicit_distance(as_polydata(mesh))
    return np.asarray(out['implicit_distance'], float)


def refine_band(sdf, origin, pitch, level, mesh, halfwidth=2.5, cap=3.0e6):
    """Replace the distance transform with EXACT distances near the isolevel.

    The transform measures to voxel CENTRES and counts a partly-covered voxel as
    solid, so it reports a radius biased outward by roughly 0.4 * pitch -- a
    sphere offset by 4.000 mm came back 4.28 mm. That bias is a constant, not
    noise, so no amount of smoothing removes it.

    Only the voxels the isosurface actually passes through matter, and that is a
    thin shell. Computing true point-to-triangle distance there costs a fraction
    of doing it everywhere and removes the bias entirely; the rest of the field
    keeps its transform value, which is only ever used to decide which side of
    the level a far-away voxel is on.
    """
    band = np.abs(sdf - level) < halfwidth * pitch
    n = int(band.sum())
    if n == 0 or n > cap:
        log('band refine: skipped (%d voxels)' % n)
        return sdf
    ii, jj, kk = np.nonzero(band)
    pts = origin + np.stack([ii, jj, kk], axis=1) * pitch
    out = sdf.copy()
    out[band] = exact_signed_distance(pts, mesh).astype(np.float32)
    log('band refine: %d voxels given exact distance' % n)
    return out


def extract(sdf, origin, pitch, level, smooth_iters):
    import pyvista as pv

    grid = pv.ImageData(dimensions=sdf.shape, spacing=(pitch, pitch, pitch),
                        origin=tuple(origin))
    grid.point_data['sdf'] = sdf.flatten(order='F')
    iso = grid.contour([float(level)], scalars='sdf')
    if iso.n_points == 0:
        raise RuntimeError('the isosurface at %.3f mm is empty -- the offset is '
                           'probably larger than the part' % level)
    iso = iso.connectivity('largest')
    iso = iso.clean().triangulate()
    if smooth_iters > 0:
        # Taubin, not Laplacian: a Laplacian smooth of a closed surface shrinks
        # it, which would silently eat into the offset we were asked for
        iso = iso.smooth_taubin(n_iter=int(smooth_iters), pass_band=0.1,
                                normalize_coordinates=True)
    return iso


def achieved(iso, mesh, sample=4000):
    """What offset we actually got, measured back against the source mesh."""
    pts = np.asarray(iso.vertices, float)
    if len(pts) > sample:
        step = len(pts) // sample
        pts = pts[::step]
    d = exact_signed_distance(pts, mesh)
    return {
        'mean_mm': float(d.mean()),
        'min_mm': float(d.min()),
        'max_mm': float(d.max()),
        'rms_err_mm': float(np.sqrt(((d - d.mean()) ** 2).mean())),
        'points': int(len(pts)),
    }


def run(job):
    import trimesh

    t0 = time.time()
    coords = np.asarray(job['coords'], dtype=float).reshape(-1, 3) * CM
    idx = np.asarray(job['idx'], dtype=np.int64).reshape(-1, 3)
    offset = float(job['offset_mm'])
    close_mm = float(job.get('close_mm', 0.0))
    smooth_iters = int(job.get('smooth_iters', 15))
    budget = float(job.get('budget', DEFAULT_BUDGET))

    mesh = trimesh.Trimesh(vertices=coords, faces=idx, process=True)
    mesh.remove_unreferenced_vertices()
    if not mesh.is_watertight:
        mesh.fill_holes()
    log('source: %d triangles, watertight %s, volume %.1f mm3'
        % (len(mesh.faces), mesh.is_watertight, mesh.volume))

    size = mesh.extents
    pitch = pick_pitch(size, offset, float(job.get('pitch_mm', 0.0) or 0.0), budget)
    margin = abs(offset) + max(close_mm, 0.0) + 4.0 * pitch
    origin = mesh.bounds[0] - margin
    span = (mesh.bounds[1] + margin) - origin
    dims = tuple(int(np.ceil(s / pitch)) + 1 for s in span)
    log('grid: pitch %.4f mm, dims %d x %d x %d (%.1f M voxels)'
        % (pitch, dims[0], dims[1], dims[2], np.prod(dims) / 1.0e6))

    solid = occupancy(mesh, origin, dims, pitch)
    log('occupancy: %.1f %% of the grid is solid' % (100.0 * solid.mean()))

    if close_mm > 0:
        solid = close_gaps(solid, pitch, close_mm)
        log('gap closing: %.3f mm applied' % close_mm)

    sdf = signed_distance(solid, pitch)
    sdf = refine_band(sdf, origin, pitch, offset, mesh)
    iso = extract(sdf, origin, pitch, offset, smooth_iters)
    log('isosurface: %d points, %d triangles' % (iso.n_points, iso.n_cells))

    faces = np.asarray(iso.faces).reshape(-1, 4)[:, 1:]
    pts = np.asarray(iso.points, float)
    out = trimesh.Trimesh(vertices=pts, faces=faces, process=True)
    out.remove_unreferenced_vertices()
    if not out.is_watertight:
        out.fill_holes()
    out.fix_normals()
    log('offset body: %d triangles, watertight %s, volume %.1f mm3'
        % (len(out.faces), out.is_watertight, out.volume))

    stats = {
        'pitch_mm': pitch,
        'dims': list(dims),
        'voxels': int(np.prod(dims)),
        'requested_offset_mm': offset,
        'close_mm': close_mm,
        'smooth_iters': smooth_iters,
        'source_tris': int(len(mesh.faces)),
        'source_volume_mm3': float(mesh.volume),
        'out_tris': int(len(out.faces)),
        'out_volume_mm3': float(out.volume),
        'watertight': bool(out.is_watertight),
        'achieved': achieved(out, mesh),
        'secs': round(time.time() - t0, 2),
    }
    a = stats['achieved']
    log('achieved offset: mean %.4f  min %.4f  max %.4f  rms err %.4f mm'
        % (a['mean_mm'], a['min_mm'], a['max_mm'], a['rms_err_mm']))

    return {
        'ok': True,
        'coords': (out.vertices * MM).ravel().tolist(),
        'idx': out.faces.ravel().astype(int).tolist(),
        'stats': stats,
        'log': list(_LOG),
    }


def main():
    jp, op = sys.argv[1], sys.argv[2]
    with open(jp, encoding='utf-8') as fh:
        job = json.load(fh)
    try:
        res = run(job)
    except Exception:
        import traceback
        res = {'ok': False, 'error': traceback.format_exc(), 'log': list(_LOG)}
    with open(op, 'w', encoding='utf-8') as fh:
        json.dump(res, fh)


if __name__ == '__main__':
    main()
