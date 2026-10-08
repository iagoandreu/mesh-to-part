"""Organic conversion engine - runs OUTSIDE Fusion, under py -3.10.

Fusion's Organic mesh conversion needs the Product Design Extension. This is a
replacement built from scratch: it wraps a cube-sphere quad cage onto the mesh,
fits one bicubic NURBS patch per cube face, and hands the six patches plus their
shared boundary curves back so the add-in can assemble a real B-Rep solid.

Why a cube-sphere and not a UV sphere: no pole singularity and no periodic seam,
so every patch is a plain open bicubic surface and the topology is 6 faces,
12 edges, 8 vertices - the simplest watertight quad layout that exists.

Why the boundaries come out watertight: adjacent patches interpolate the *same*
row of wrapped points, with the same uniform parameterisation and the same end
conditions, so their shared boundary curve is identical by construction rather
than by tolerance.

Fusion's bundled interpreter has no numpy, hence the split. All lengths are
centimetres, matching Fusion's internal unit.

    py -3.10 mtp_organic.py job.json out.json
"""

import json
import math
import sys
import time

import numpy as np
from scipy.interpolate import BSpline
from scipy.spatial import cKDTree

import trimesh


# --------------------------------------------------------------------------- #
# exact closest point on a triangle soup
#
# trimesh.proximity needs rtree, which is not installed here and is not worth
# adding to the user's interpreter for this. A KD-tree over triangle centroids
# shortlists candidates and the exact point-triangle projection finishes the job;
# the shortlist radius is bounded so the answer is exact, not approximate.
# --------------------------------------------------------------------------- #

class Projector(object):
    def __init__(self, mesh):
        self.tri = np.asarray(mesh.triangles, float)          # (F, 3, 3)
        self.a = self.tri[:, 0, :]
        self.b = self.tri[:, 1, :]
        self.c = self.tri[:, 2, :]
        self.centroid = self.tri.mean(axis=1)
        self.radius = float(np.linalg.norm(
            self.tri - self.centroid[:, None, :], axis=2).max())
        self.tree = cKDTree(self.centroid)
        n = np.cross(self.b - self.a, self.c - self.a)
        ln = np.linalg.norm(n, axis=1)
        self.normal = n / np.maximum(ln, 1e-30)[:, None]

    def __call__(self, points, k=24, exact=True):
        points = np.asarray(points, float)
        k = min(k, len(self.centroid))
        d0, idx0 = self.tree.query(points, k=k, workers=-1)
        if idx0.ndim == 1:
            idx0 = idx0[:, None]
        out_pts, out_d, out_f = self._batch(points, idx0)
        # anything whose centroid is farther than this bound cannot win
        bound = out_d + self.radius
        if not exact:
            return out_pts, out_d, out_f
        far = d0[:, -1] < bound if d0.ndim > 1 else np.zeros(len(points), bool)
        for i in np.nonzero(far)[0]:
            cand = self.tree.query_ball_point(points[i], bound[i], workers=-1)
            if len(cand) <= k:
                continue
            p, d, f = self._batch(points[i][None, :], np.array(cand)[None, :])
            if d[0] < out_d[i]:
                out_pts[i] = p[0]
                out_d[i] = d[0]
                out_f[i] = f[0]
        return out_pts, out_d, out_f

    def _batch(self, points, cand):
        n, m = cand.shape
        p = np.repeat(points, m, axis=0)
        f = cand.reshape(-1)
        q = closest_on_triangle(p, self.a[f], self.b[f], self.c[f])
        d = np.linalg.norm(q - p, axis=1).reshape(n, m)
        j = d.argmin(axis=1)
        rows = np.arange(n)
        return q.reshape(n, m, 3)[rows, j], d[rows, j], cand[rows, j]


def closest_on_triangle(p, a, b, c):
    """Vectorised closest point on a triangle (Ericson, Real-Time Collision
    Detection 5.1.5). Regions are assigned in reverse priority so the highest
    priority case is the one that survives."""
    ab = b - a
    ac = c - a
    ap = p - a
    d1 = (ab * ap).sum(1)
    d2 = (ac * ap).sum(1)
    bp = p - b
    d3 = (ab * bp).sum(1)
    d4 = (ac * bp).sum(1)
    cp = p - c
    d5 = (ab * cp).sum(1)
    d6 = (ac * cp).sum(1)
    va = d3 * d6 - d5 * d4
    vb = d5 * d2 - d1 * d6
    vc = d1 * d4 - d3 * d2
    denom = va + vb + vc
    safe = np.where(np.abs(denom) < 1e-20, 1.0, denom)

    v = vb / safe
    w = vc / safe
    out = a + ab * v[:, None] + ac * w[:, None]              # interior

    # edge BC
    d43 = d4 - d3
    d56 = d5 - d6
    t = np.where(np.abs(d43 - d56) < 1e-20, 0.0, d43 / np.where(
        np.abs(d43 - d56) < 1e-20, 1.0, d43 - d56))
    t = np.clip(t, 0.0, 1.0)
    mBC = (va <= 0) & (d43 >= 0) & (d56 >= 0)
    out = np.where(mBC[:, None], b + (c - b) * t[:, None], out)

    # edge AC
    tw = np.where(np.abs(d2 - d6) < 1e-20, 0.0,
                  d2 / np.where(np.abs(d2 - d6) < 1e-20, 1.0, d2 - d6))
    tw = np.clip(tw, 0.0, 1.0)
    mAC = (vb <= 0) & (d2 >= 0) & (d6 <= 0)
    out = np.where(mAC[:, None], a + ac * tw[:, None], out)

    # edge AB
    tv = np.where(np.abs(d1 - d3) < 1e-20, 0.0,
                  d1 / np.where(np.abs(d1 - d3) < 1e-20, 1.0, d1 - d3))
    tv = np.clip(tv, 0.0, 1.0)
    mAB = (vc <= 0) & (d1 >= 0) & (d3 <= 0)
    out = np.where(mAB[:, None], a + ab * tv[:, None], out)

    out = np.where(((d6 >= 0) & (d5 <= d6))[:, None], c, out)
    out = np.where(((d3 >= 0) & (d4 <= d3))[:, None], b, out)
    out = np.where(((d1 <= 0) & (d2 <= 0))[:, None], a, out)
    return out


# --------------------------------------------------------------------------- #
# the cage
# --------------------------------------------------------------------------- #

# Each cube face as (u_dir, v_dir, normal) with u x v == outward normal, so a
# loop running u+ then v+ is counter-clockwise seen from outside and the patch
# normal needs no flipping later.
CUBE_FACES = [
    ((0, 1, 0), (0, 0, 1), (1, 0, 0)),      # +X
    ((0, 0, 1), (0, 1, 0), (-1, 0, 0)),     # -X
    ((0, 0, 1), (1, 0, 0), (0, 1, 0)),      # +Y
    ((1, 0, 0), (0, 0, 1), (0, -1, 0)),     # -Y
    ((1, 0, 0), (0, 1, 0), (0, 0, 1)),      # +Z
    ((0, 1, 0), (1, 0, 0), (0, 0, -1)),     # -Z
]


def cube_sphere(n):
    """Cube-sphere cage with n+1 points per side per face.

    Returns (verts, grids) where grids[f] is an (n+1, n+1) array of indices into
    verts. Points on shared cube edges and corners resolve to the SAME index, so
    the cage is one welded surface and the wrap can never tear it open.
    """
    t = np.linspace(-1.0, 1.0, n + 1)
    lut = {}
    verts = []
    grids = []
    for u_dir, v_dir, normal in CUBE_FACES:
        u_dir = np.array(u_dir, float)
        v_dir = np.array(v_dir, float)
        normal = np.array(normal, float)
        g = np.zeros((n + 1, n + 1), dtype=np.int64)
        for i in range(n + 1):
            for j in range(n + 1):
                p = normal + u_dir * t[i] + v_dir * t[j]
                key = (round(p[0], 9), round(p[1], 9), round(p[2], 9))
                idx = lut.get(key)
                if idx is None:
                    idx = len(verts)
                    lut[key] = idx
                    verts.append(p)
                g[i, j] = idx
        grids.append(g)
    verts = np.array(verts, float)
    verts /= np.linalg.norm(verts, axis=1)[:, None]     # cube -> sphere
    return verts, grids


def cage_neighbours(n_verts, grids):
    """Neighbour lists from the quad edges, used for Laplacian smoothing."""
    nb = [set() for _ in range(n_verts)]
    for g in grids:
        n = g.shape[0]
        for i in range(n):
            for j in range(n):
                a = g[i, j]
                if i + 1 < n:
                    b = g[i + 1, j]
                    nb[a].add(b)
                    nb[b].add(a)
                if j + 1 < n:
                    b = g[i, j + 1]
                    nb[a].add(b)
                    nb[b].add(a)
    width = max(len(s) for s in nb)
    table = np.full((n_verts, width), -1, dtype=np.int64)
    for i, s in enumerate(nb):
        for k, v in enumerate(sorted(s)):
            table[i, k] = v
    return table


def laplacian(points, table):
    valid = table >= 0
    safe = np.where(valid, table, 0)
    acc = (points[safe] * valid[:, :, None]).sum(axis=1)
    cnt = valid.sum(axis=1)[:, None].astype(float)
    return acc / np.maximum(cnt, 1.0) - points


def smooth(points, table, weight, normals=None):
    """One Laplacian pass. With normals supplied only the TANGENTIAL part of the
    displacement is applied, which is what lets the cage settle into concavities:
    plain Laplacian smoothing lifts points off the surface and the projection then
    drags them back to the nearest convex bit, so a hollow never gets filled."""
    if weight <= 0.0:
        return points
    disp = laplacian(points, table)
    if normals is not None:
        disp = disp - (disp * normals).sum(axis=1)[:, None] * normals
    return points + weight * disp


# --------------------------------------------------------------------------- #
# the wrap
# --------------------------------------------------------------------------- #

def neighbour_mean(values, table):
    valid = table >= 0
    safe = np.where(valid, table, 0)
    acc = (values[safe] * valid[:, :, None]).sum(axis=1)
    cnt = valid.sum(axis=1)[:, None].astype(float)
    return acc / np.maximum(cnt, 1.0)


def spacing(points, table):
    """Mean distance from each point to its cage neighbours."""
    valid = table >= 0
    safe = np.where(valid, table, 0)
    d = np.linalg.norm(points[safe] - points[:, None, :], axis=2) * valid
    return d.sum(axis=1) / np.maximum(valid.sum(axis=1), 1)


def spacing_cv(points, table):
    """Coefficient of variation of cage edge length. This, not the Laplacian
    magnitude, is what uniform-knot fitting actually cares about: a curved cage
    has a large Laplacian even when it is perfectly evenly spaced."""
    valid = table >= 0
    safe = np.where(valid, table, 0)
    d = np.linalg.norm(points[safe] - points[:, None, :], axis=2)[valid]
    return float(d.std() / max(d.mean(), 1e-30))


def quad_arrays(grids):
    """Flat arrays of the four corner indices of every cage quad, once."""
    a, b, c, d = [], [], [], []
    for g in grids:
        n = g.shape[0] - 1
        a.append(g[:n, :n].ravel())
        b.append(g[1:, :n].ravel())      # +u
        c.append(g[1:, 1:].ravel())      # +u +v
        d.append(g[:n, 1:].ravel())      # +v
    return (np.concatenate(a), np.concatenate(b),
            np.concatenate(c), np.concatenate(d))


def folded_quads(points, quads, proj=None):
    """Which cage quads are turned inside out relative to their own neighbourhood.

    Comparing a quad against the MESH normal underneath it looked like the obvious
    test and it is wrong on a creased scan: a quad straddling a 150-degree crease
    projects its centre onto one side or the other, and the mesh normal there can
    oppose the quad's own average direction with nothing actually folded. That test
    reported 2% of quads folded on a body whose volume was right to 0.16% and which
    healed into a valid solid - false positives.

    The intrinsic test compares each quad with the area-weighted vertex normals of
    its own corners. A genuinely buckled quad opposes the field its neighbours
    define; a quad sitting honestly across a crease does not. It needs no projection
    either, so it is cheap enough to run inside every step.
    """
    a, b, c, d = quads
    nq = np.cross(points[b] - points[a], points[d] - points[a])
    ln = np.linalg.norm(nq, axis=1)
    nq_unit = nq / np.maximum(ln, 1e-30)[:, None]

    vn = np.zeros_like(points)
    for arr in (a, b, c, d):
        np.add.at(vn, arr, nq)              # area weighted: nq is not normalised
    field = vn[a] + vn[b] + vn[c] + vn[d]
    fl = np.linalg.norm(field, axis=1)
    field = field / np.maximum(fl, 1e-30)[:, None]
    return (nq_unit * field).sum(axis=1) < 0.0, ln


def quad_metrics(points, quads):
    """Per-quad aspect ratio and skew.

    quad_arrays orders the corners a -> b (+u) -> c (+u+v) -> d (+v), so the two
    u-edges are a-b and d-c, and the two v-edges are a-d and b-c.

    aspect  longer edge family over the shorter; 1.0 is square
    skew    |cos| of the angle between the u and v directions; 0.0 is square,
            1.0 is a fully collapsed sliver
    """
    a, b, c, d = quads
    ab = points[b] - points[a]
    dc = points[c] - points[d]
    ad = points[d] - points[a]
    bc = points[c] - points[b]
    eu = 0.5 * (np.linalg.norm(ab, axis=1) + np.linalg.norm(dc, axis=1))
    ev = 0.5 * (np.linalg.norm(ad, axis=1) + np.linalg.norm(bc, axis=1))
    aspect = np.maximum(eu, ev) / np.maximum(np.minimum(eu, ev), 1.0e-30)
    u = ab + dc
    v = ad + bc
    u = u / np.maximum(np.linalg.norm(u, axis=1), 1.0e-30)[:, None]
    v = v / np.maximum(np.linalg.norm(v, axis=1), 1.0e-30)[:, None]
    skew = np.abs((u * v).sum(axis=1))
    return aspect, skew, eu, ev


def aspect_report(points, quads, limit):
    aspect, skew, eu, ev = quad_metrics(points, quads)
    L = edge_lengths(points, quads)
    ratio = L.max(axis=1) / np.maximum(L.min(axis=1), 1.0e-30)
    return {
        # the limit applies to the EDGE RATIO within a quad -- that is the
        # property that decides whether an element is usable
        'ratio_median': float(np.median(ratio)),
        'ratio_p95': float(np.percentile(ratio, 95)),
        'ratio_max': float(ratio.max()),
        'over_limit': int((ratio > limit).sum()),
        'min_edge_mm': float(L.min()) * 10.0,
        'median_edge_mm': float(np.median(L)) * 10.0,
        'aspect_median': float(np.median(aspect)),
        'aspect_p95': float(np.percentile(aspect, 95)),
        'aspect_max': float(aspect.max()),
        'skew_p95': float(np.percentile(skew, 95)),
        'skew_max': float(skew.max()),
        'quads': int(len(aspect)),
        'limit': float(limit),
    }


def vertex_badness(points, quads, n_verts, limit):
    """How much each cage vertex is implicated in a badly shaped quad.

    Only quads that actually breach the limit contribute anything, so a
    relaxation driven by this leaves the square, well-fitting regions of the
    cage alone rather than smoothing everything and washing out the detail that
    was already landing correctly.
    """
    ratio = edge_ratio(points, quads)
    bad = np.maximum(ratio / max(limit, 1.0e-9) - 1.0, 0.0)
    acc = np.zeros(n_verts)
    cnt = np.zeros(n_verts)
    for arr in quads:
        np.add.at(acc, arr, bad)
        np.add.at(cnt, arr, 1.0)
    return acc / np.maximum(cnt, 1.0)


def relax_quality(points, table, quads, proj, fid, iters, limit, gain, log):
    """Tangential relaxation weighted by how badly shaped each quad is.

    A uniform relaxation has to stay weak or it washes out the parts of the cage
    that are tracking the surface properly -- which leaves the long thin slivers
    that collect along a crease exactly as they were. Those slivers are what the
    limit surface renders as a crease of its own, and what makes a downstream
    fillet fail. Weighting each vertex's step by the badness of the quads it
    belongs to spends the smoothing only where the shape is actually wrong.
    """
    for it in range(iters):
        nrm = proj.normal[fid]
        nrm = nrm + neighbour_mean(nrm, table)
        nrm /= np.maximum(np.linalg.norm(nrm, axis=1), 1e-30)[:, None]

        bad = vertex_badness(points, quads, len(points), limit)
        if bad.max() <= 0.0:
            if log:
                log('  quality %3d/%d  every quad already within aspect %.2f'
                    % (it + 1, iters, limit))
            break
        w = np.clip(gain * bad, 0.0, 0.9)

        disp = laplacian(points, table)
        disp = disp - (disp * nrm).sum(axis=1)[:, None] * nrm
        step = w[:, None] * disp
        # Cap the step against the TYPICAL cage spacing, not the point's own.
        # A sliver's local spacing is tiny by definition, so a proportional cap
        # freezes precisely the points that need to travel to escape it -- the
        # first version of this pass moved the worst aspect on a real scan from
        # 47.7 only to 35.7 for exactly that reason.
        sp = spacing(points, table)
        lim = 0.35 * np.maximum(sp, np.median(sp))
        ln = np.linalg.norm(step, axis=1)
        step = step * np.minimum(1.0, lim / np.maximum(ln, 1e-30))[:, None]
        points, _bad = safe_step_mesh(points, step, quads, proj)
        points, _d, fid = proj(points, k=16, exact=False)
        if log and (it % 10 == 0 or it == iters - 1):
            r = aspect_report(points, quads, limit)
            log('  quality %3d/%d  edge ratio med %.2f p95 %.2f max %.2f'
                '  over limit %d'
                % (it + 1, iters, r['ratio_median'], r['ratio_p95'],
                   r['ratio_max'], r['over_limit']))
    return points, fid


def square_targets(points, quads):
    """Where each quad's four corners would sit if the quad were a square.

    The aspect-weighted Laplacian attacks unequal EDGE LENGTHS only. It left
    skew at p95 0.94 on a real scan, because a Laplacian has no notion of a
    corner angle: a perfectly parallelogram-shaped sliver has a small Laplacian
    and is still a sliver. This target attacks both at once.

    The square is built in the quad's own best-fit plane, at its own centroid
    and mean corner radius, rotated to whichever of the four 90-degree
    orientations best matches the quad as it stands (a circular mean of the
    de-rotated corner angles). So the target is the NEAREST square, which keeps
    the cage where it is instead of spinning it.
    """
    a, b, c, d = quads
    P = (points[a], points[b], points[c], points[d])
    g = (P[0] + P[1] + P[2] + P[3]) * 0.25
    n = np.cross(P[2] - P[0], P[3] - P[1])
    n = n / np.maximum(np.linalg.norm(n, axis=1), 1.0e-30)[:, None]
    e1 = P[0] - g
    e1 = e1 - (e1 * n).sum(axis=1)[:, None] * n
    e1 = e1 / np.maximum(np.linalg.norm(e1, axis=1), 1.0e-30)[:, None]
    e2 = np.cross(n, e1)

    r = np.zeros(len(g))
    cs = np.zeros(len(g))
    sn = np.zeros(len(g))
    for k in range(4):
        v = P[k] - g
        x = (v * e1).sum(axis=1)
        y = (v * e2).sum(axis=1)
        r += np.sqrt(x * x + y * y)
        t = np.arctan2(y, x) - k * (np.pi * 0.5)
        cs += np.cos(t)
        sn += np.sin(t)
    r *= 0.25
    th0 = np.arctan2(sn, cs)

    out = []
    for k in range(4):
        t = th0 + k * (np.pi * 0.5)
        out.append(g + (r * np.cos(t))[:, None] * e1
                   + (r * np.sin(t))[:, None] * e2)
    return out


def square_displacement(points, quads):
    """Per-vertex move toward the average of its incident quads' square targets."""
    tgt = square_targets(points, quads)
    acc = np.zeros_like(points)
    cnt = np.zeros(len(points))
    for arr, t in zip(quads, tgt):
        np.add.at(acc, arr, t)
        np.add.at(cnt, arr, 1.0)
    return acc / np.maximum(cnt, 1.0)[:, None] - points


def min_edge_push(points, quads, min_edge):
    """Separate the endpoints of any cage edge shorter than min_edge.

    Squareness alone permits a quad that is square and tiny. A Form element
    below a millimetre or so is unusable for filleting whatever its shape, so
    the floor is enforced directly rather than hoped for.
    """
    a, b, c, d = quads
    acc = np.zeros_like(points)
    cnt = np.zeros(len(points))
    for i, j in ((a, b), (b, c), (c, d), (d, a)):
        e = points[j] - points[i]
        L = np.linalg.norm(e, axis=1)
        short = L < min_edge
        if not short.any():
            continue
        u = e / np.maximum(L, 1.0e-30)[:, None]
        push = (np.where(short, 0.5 * (min_edge - L), 0.0))[:, None] * u
        np.add.at(acc, i, -push)
        np.add.at(cnt, i, 1.0)
        np.add.at(acc, j, push)
        np.add.at(cnt, j, 1.0)
    return acc / np.maximum(cnt, 1.0)[:, None]


def valence_damping(table, rings=2):
    """1 in the regular interior, 0 at the cage's extraordinary vertices.

    A cube-sphere cage has 8 valence-3 corners. "Nearest square" is not defined
    where only three quads meet, and pushing toward it there is what Fusion
    rejects as ASM_TSP_IG_SURFACE_PARALLEL_SINGULARITY -- the T-Spline
    validator sees degenerate, near-parallel edges at the singular point. The
    ramp keeps the force continuous instead of switching it off at a hard
    boundary, which would just move the problem one ring out.
    """
    valid = table >= 0
    safe = np.where(valid, table, 0)
    cnt = np.maximum(valid.sum(axis=1), 1)
    damp = (valid.sum(axis=1) == 4).astype(float)
    for _ in range(int(rings)):
        nb = (damp[safe] * valid).sum(axis=1) / cnt
        damp = np.minimum(damp, nb)
    return damp


def edge_lengths(points, quads):
    """The four edge lengths of every cage quad, in order a-b, b-c, c-d, d-a."""
    a, b, c, d = quads
    return np.stack([
        np.linalg.norm(points[b] - points[a], axis=1),
        np.linalg.norm(points[c] - points[b], axis=1),
        np.linalg.norm(points[d] - points[c], axis=1),
        np.linalg.norm(points[a] - points[d], axis=1)], axis=1)


def edge_ratio(points, quads):
    """Longest edge over shortest edge WITHIN each quad.

    This, not the aspect of the two edge families and not the corner angle, is
    what decides whether an element is usable: an element with one edge five
    times another is the thing that shows as a crease and defeats a fillet. A
    quad whose corners are not 90 degrees is fine as long as its four edges are
    close to equal.
    """
    L = edge_lengths(points, quads)
    return L.max(axis=1) / np.maximum(L.min(axis=1), 1.0e-30)


def edge_equalise(points, quads):
    """Spring force pulling the four edges of each quad toward a common length.

    Deliberately NOT the nearest-square target that this replaced. That one also
    rotated quads toward 90-degree corners, and the rotation is what Fusion
    refused as ASM_TSP_IG_SURFACE_PARALLEL_SINGULARITY -- it moves the cage
    across the surface. Equalising lengths only redistributes points along the
    grid lines they already lie on, which is a far smaller perturbation for the
    same benefit to element usability.
    """
    a, b, c, d = quads
    L = edge_lengths(points, quads)
    rest = L.mean(axis=1)
    acc = np.zeros_like(points)
    cnt = np.zeros(len(points))
    for k, (i, j) in enumerate(((a, b), (b, c), (c, d), (d, a))):
        e = points[j] - points[i]
        ln = np.maximum(np.linalg.norm(e, axis=1), 1.0e-30)
        u = e / ln[:, None]
        # positive when this edge is longer than its quad's mean: the endpoints
        # move toward each other, and a short edge pushes them apart
        move = (0.5 * (ln - rest))[:, None] * u
        np.add.at(acc, i, move)
        np.add.at(cnt, i, 1.0)
        np.add.at(acc, j, -move)
        np.add.at(cnt, j, 1.0)
    return acc / np.maximum(cnt, 1.0)[:, None]


def relax_square(points, table, quads, proj, fid, iters, weight, min_edge,
                 limit, log, step_cap=0.10):
    damp = valence_damping(table, rings=1)
    """Drive the cage toward square elements, tangentially, with an edge floor.

    Tangential so the cage stays on the mesh: the whole point of the wrap is
    that the cage sits on the surface, and a squareness force with a normal
    component would lift it off and the projection would drag it back to the
    nearest convex face, undoing the concavity work.
    """
    for it in range(iters):
        nrm = proj.normal[fid]
        nrm = nrm + neighbour_mean(nrm, table)
        nrm /= np.maximum(np.linalg.norm(nrm, axis=1), 1e-30)[:, None]

        disp = weight * edge_equalise(points, quads)
        if min_edge > 0.0:
            disp = disp + min_edge_push(points, quads, min_edge)
        disp = disp * damp[:, None]
        disp = disp - (disp * nrm).sum(axis=1)[:, None] * nrm

        # A shape pass must NUDGE, not transport. At 0.35 * spacing per pass,
        # 100 passes let a point travel far enough to slide around a toe and
        # put the cage through itself -- Fusion then refuses the T-Spline with
        # ASM_TSP_IG_SURFACE_PARALLEL_SINGULARITY, which the intrinsic fold test
        # cannot see because the two crossing quads are not neighbours.
        sp = spacing(points, table)
        lim = step_cap * float(np.median(sp))
        ln = np.linalg.norm(disp, axis=1)
        disp = disp * np.minimum(1.0, lim / np.maximum(ln, 1e-30))[:, None]
        points, _bad = safe_step_mesh(points, disp, quads, proj)
        points, _d, fid = proj(points, k=16, exact=False)
        if log and (it % 20 == 0 or it == iters - 1):
            r = aspect_report(points, quads, limit)
            log('  edges %3d/%d  ratio med %.2f p95 %.2f max %.2f  over %d'
                '  shortest %.3f mm'
                % (it + 1, iters, r['ratio_median'], r['ratio_p95'],
                   r['ratio_max'], r['over_limit'], r['min_edge_mm']))
    return points, fid


def safe_step_mesh(points, step, quads, proj, tries=4, bad_before=None):
    """Apply a displacement, backing off only where it would create a NEW fold.

    Two things had to be right here. The cheap self-consistency test (did this
    quad's normal flip during this step?) misses buckles that build up gradually,
    which is exactly how they formed - the count climbed steadily while no single
    step ever flipped anything. And backing off around quads that were ALREADY
    folded starves those vertices of movement: with that bug the cage finished
    72 mm short of the surface in places. So the comparison is against the surface
    normal underneath, and only against folds this step is responsible for.
    """
    a, b, c, d = quads
    if bad_before is None:
        bad_before, _l = folded_quads(points, quads, proj)
    scale = np.ones(len(points))
    cand = points + step
    for _ in range(tries):
        bad, _ln = folded_quads(cand, quads, proj)
        fresh = bad & ~bad_before
        if not fresh.any():
            return cand, bad
        touched = np.zeros(len(points), bool)
        for arr in (a, b, c, d):
            touched[arr[fresh]] = True
        scale[touched] *= 0.4
        cand = points + step * scale[:, None]
    bad, _ln = folded_quads(cand, quads, proj)
    return cand, bad


def untangle(points, table, quads, proj, passes, log):
    """Iteratively unbuckle the cage, giving up accuracy locally to do it.

    A folded quad is fatal for the Form output: Fusion rejects a self-intersecting
    T-Spline with ASM_TSP_IG_SURFACE_PARALLEL_SINGULARITY and produces no body at
    all. So this is not a quality nicety, it is a hard requirement.

    The earlier version smoothed the folded region and then re-projected it onto
    the mesh, which is why it stalled at 65 folds: the projection put the buckle
    straight back. Here the folded neighbourhood is smoothed and deliberately left
    OFF the surface, while everything else stays projected onto it. Trading a
    little local fidelity for a body that exists is the right way round.
    """
    a, b, c, d = quads
    n_before = None
    weight = 0.5
    for it in range(passes):
        bad, _ln = folded_quads(points, quads, proj)
        n_bad = int(bad.sum())
        if n_before is None:
            n_before = n_bad
        if n_bad == 0:
            break
        touched = np.zeros(len(points), bool)
        for arr in (a, b, c, d):
            touched[arr[bad]] = True
        for _ring in range(2):                 # widen so the repair is not a dent
            nb = table[touched]
            touched[nb[nb >= 0]] = True

        disp = laplacian(points, table)
        points = points + (weight * touched[:, None]) * disp
        free = ~touched
        if free.any():
            pf, _dd, _ff = proj(points[free], k=16, exact=False)
            points[free] = pf
        weight = min(0.9, weight * 1.15)       # push harder if it is not clearing
        if log and (it % 5 == 0 or it == passes - 1):
            log('  untangle %2d/%d  folded quads %d  weight %.2f'
                % (it + 1, passes, n_bad, weight))
    bad, _ln = folded_quads(points, quads, proj)
    if log:
        log('  untangle: folded quads %s -> %d' % (n_before, int(bad.sum())))
    return points, int(bad.sum())


def shrink_wrap(mesh, proj, points, table, grids, iters, smooth_start, smooth_end,
                sample_count, tangential_iters, untangle_passes, log,
                aspect_limit=0.0, aspect_iters=0, aspect_gain=1.2,
                square_iters=0, square_weight=0.5, min_edge=0.0,
                square_step_cap=0.10):
    """Pull the cage onto the mesh, then let it slide along the surface.

    Phase 1 - free wrap. A KD-tree over a dense surface sample drives it, which is
    far cheaper than exact projection, and full Laplacian smoothing keeps the cage
    from folding while it travels in from the enclosing sphere.

    Phase 2 - tangential relaxation. Every point sits exactly on the mesh and only
    slides along it. This is the phase that reaches into concavities: normal-
    direction smoothing would lift a point out of a hollow and the projection would
    then snap it back to the nearest convex face, so the hollow could never fill.
    """
    samples, _ = trimesh.sample.sample_surface(mesh, sample_count)
    tree = cKDTree(samples)
    quads = quad_arrays(grids)

    for it in range(iters):
        frac = it / max(1, iters - 1)
        w = smooth_start + (smooth_end - smooth_start) * frac
        _, idx = tree.query(points, k=1, workers=-1)
        pull = 0.85 * (samples[idx] - points)
        # A cage row that overshoots into a concavity crosses its neighbours and
        # the cage buckles; most folds were being created right here, before the
        # relaxation ever ran. Cap the advance at the local spacing and refuse any
        # step that inverts a quad.
        lim = 1.0 * spacing(points, table)
        ln = np.linalg.norm(pull, axis=1)
        pull *= np.minimum(1.0, lim / np.maximum(ln, 1e-30))[:, None]
        points, bad_now = safe_step_mesh(points, pull, quads, proj)
        sm = smooth(points, table, w) - points
        points, bad_now = safe_step_mesh(points, sm, quads, proj,
                                         bad_before=bad_now)
        if log and (it % 10 == 0 or it == iters - 1):
            d = np.linalg.norm(samples[idx] - points, axis=1)
            bad, _l = folded_quads(points, quads, proj)
            log('  wrap %3d/%d  smooth %.3f  residual max %.3f mean %.3f mm  folded %d'
                % (it + 1, iters, w, d.max() * 10.0, d.mean() * 10.0, int(bad.sum())))

    points, _d, fid = proj(points, k=16, exact=False)
    bad0, _l = folded_quads(points, quads, proj)
    if log:
        log('  folded quads after the wrap: %d of %d' % (int(bad0.sum()), len(quads[0])))
    for it in range(tangential_iters):
        w = smooth_end + (smooth_start - smooth_end) *             (1.0 - it / max(1, tangential_iters - 1)) * 0.5
        # Smooth the surface normal over the cage neighbourhood first. The raw
        # face normal flips across a crease, which sends the tangential step in a
        # different direction on every pass and makes the relaxation oscillate
        # instead of settle.
        nrm = proj.normal[fid]
        nrm = nrm + neighbour_mean(nrm, table)
        nrm /= np.maximum(np.linalg.norm(nrm, axis=1), 1e-30)[:, None]

        disp = laplacian(points, table)
        disp = disp - (disp * nrm).sum(axis=1)[:, None] * nrm
        step = w * disp
        # never move a point more than a fraction of its own spacing in one pass
        lim = 0.35 * spacing(points, table)
        ln = np.linalg.norm(step, axis=1)
        scale = np.minimum(1.0, lim / np.maximum(ln, 1e-30))
        points, _bad = safe_step_mesh(points, step * scale[:, None], quads, proj)
        points, _d, fid = proj(points, k=16, exact=False)
        if log and (it % 15 == 0 or it == tangential_iters - 1):
            bad, _l = folded_quads(points, quads, proj)
            log('  relax %3d/%d  w %.3f  spacing CV %.4f  folded %d'
                % (it + 1, tangential_iters, w, spacing_cv(points, table),
                   int(bad.sum())))
    quality = None
    if aspect_iters > 0 and aspect_limit > 0.0:
        before = aspect_report(points, quads, aspect_limit)
        if log:
            log('  cage quads before quality relax: edge ratio med %.2f '
                'p95 %.2f max %.2f, %d of %d over %.2f'
                % (before['ratio_median'], before['ratio_p95'],
                   before['ratio_max'], before['over_limit'],
                   before['quads'], aspect_limit))
        points, fid = relax_quality(points, table, quads, proj, fid,
                                    int(aspect_iters), float(aspect_limit),
                                    float(aspect_gain), log)
        quality = {'before': before,
                   'after': aspect_report(points, quads, aspect_limit)}

    if square_iters > 0:
        lim_for_report = aspect_limit if aspect_limit > 0 else 3.0
        if quality is None:
            quality = {'before': aspect_report(points, quads, lim_for_report)}
        points, fid = relax_square(points, table, quads, proj, fid,
                                   int(square_iters), float(square_weight),
                                   float(min_edge), lim_for_report, log,
                                   step_cap=float(square_step_cap))

    points, n_folded = untangle(points, table, quads, proj,
                                int(untangle_passes), log)
    points, _d, _f = proj(points)          # exact landing
    rep_limit = aspect_limit if aspect_limit > 0.0 else 3.0
    if quality is not None:
        # untangle moves points too, so the number that gets reported is the one
        # measured on the cage that is actually used
        quality['after'] = aspect_report(points, quads, rep_limit)
        _a, _s, eu, ev = quad_metrics(points, quads)
        quality['after']['min_edge_mm'] = float(min(eu.min(), ev.min())) * 10.0
    elif aspect_limit > 0.0:
        quality = {'after': aspect_report(points, quads, rep_limit)}
    return points, n_folded, quality


# --------------------------------------------------------------------------- #
# the patches
# --------------------------------------------------------------------------- #

def clamped_knots(count, degree=3):
    """Clamped uniform knot vector for `count` control points."""
    inner = count - degree - 1
    if inner < 0:
        raise ValueError('need at least %d control points' % (degree + 1))
    return np.concatenate([
        np.zeros(degree + 1),
        np.arange(1, inner + 1, dtype=float) / (inner + 1),
        np.ones(degree + 1),
    ])


def design(count, m, degree=3):
    """Basis matrix (m x count) at m uniform parameters, plus the knot vector."""
    kv = clamped_knots(count, degree)
    x = np.linspace(0.0, 1.0, m)
    x[-1] = np.nextafter(1.0, 0.0)          # keep the last point inside the span
    B = BSpline.design_matrix(x, kv, degree).toarray()
    return B, kv


def second_difference(count):
    """Fairness operator. Penalising the second difference of the control net is
    what stops a least-squares fit from ringing between samples."""
    if count < 3:
        return np.zeros((0, count))
    D = np.zeros((count - 2, count))
    for i in range(count - 2):
        D[i, i] = 1.0
        D[i, i + 1] = -2.0
        D[i, i + 2] = 1.0
    return D


def fit_curve_fixed_ends(data, B, lam):
    """Least-squares cubic B-spline through a row of points, with the two end
    control points pinned to the end data points.

    Pinning the ends is what makes three patches meeting at a cage corner agree on
    that corner exactly, and computing this from the shared data row alone is what
    makes two patches agree on the whole shared boundary.
    """
    m, P = B.shape
    c = np.zeros((P, 3))
    c[0] = data[0]
    c[-1] = data[-1]
    free = np.arange(1, P - 1)
    if len(free) == 0:
        return c
    A = B[:, free]
    rhs = data - np.outer(B[:, 0], c[0]) - np.outer(B[:, -1], c[-1])
    D = second_difference(P)
    Df = D[:, free]
    Dk = D[:, 0:1] @ c[0][None, :] + D[:, -1:] @ c[-1][None, :]
    # Stacked least squares, not the normal equations: forming A.T@A squares the
    # condition number, and at high control density that alone was enough to put
    # metres of ringing into the fit.
    lhs = np.vstack([A, math.sqrt(lam) * Df])
    r = np.vstack([rhs, -math.sqrt(lam) * Dk])
    c[free] = np.linalg.lstsq(lhs, r, rcond=None)[0]
    return c


def fit_patch_constrained(data, B, lam, bnd_u0, bnd_u1, bnd_v0, bnd_v1):
    """Least-squares bicubic patch with all four boundary control rows given.

    data is (m, m, 3) indexed [u][v]; the returned net is (P, P, 3) indexed [u][v].
    Interior control points solve  min ||Bu C Bv^T - D||^2 + lam * fairness  with
    the boundary entries held at the values the shared edge fits produced.
    """
    m, P = B.shape
    C = np.zeros((P, P, 3))
    C[0, :, :] = bnd_u0
    C[-1, :, :] = bnd_u1
    C[:, 0, :] = bnd_v0
    C[:, -1, :] = bnd_v1

    mask = np.ones((P, P), bool)
    mask[0, :] = mask[-1, :] = mask[:, 0] = mask[:, -1] = False
    free = np.flatnonzero(mask.ravel())
    if len(free) == 0:
        return C

    K = np.kron(B, B)                            # (m*m, P*P), row = u*m+v
    D2 = second_difference(P)
    I = np.eye(P)
    S = np.vstack([np.kron(D2, I), np.kron(I, D2)])

    A = K[:, free]
    fixed = np.flatnonzero(~mask.ravel())
    known = C.reshape(-1, 3)[fixed]
    rhs = data.reshape(-1, 3) - K[:, fixed] @ known
    Sf = S[:, free]
    Sk = S[:, fixed] @ known
    w = math.sqrt(lam)
    sol = np.linalg.lstsq(np.vstack([A, w * Sf]),
                          np.vstack([rhs, -w * Sk]), rcond=None)[0]
    flat = C.reshape(-1, 3)
    flat[free] = sol
    return flat.reshape(P, P, 3)


def eval_surface(net, knots_u, knots_v, nu, nv):
    """Sample the fitted surface, so its accuracy can be measured before Fusion
    ever sees it. net is indexed [u][v]; the result is indexed [u][v] too."""
    qu = np.linspace(0.0, 1.0, nu)
    qv = np.linspace(0.0, 1.0, nv)
    qu[-1] = qv[-1] = np.nextafter(1.0, 0.0)
    Bu = BSpline.design_matrix(qu, knots_u, 3).toarray()
    Bv = BSpline.design_matrix(qv, knots_v, 3).toarray()
    out = np.empty((nu, nv, 3))
    for d in range(3):
        out[:, :, d] = Bu @ net[:, :, d] @ Bv.T
    return out


# --------------------------------------------------------------------------- #
# TSM writer - the real Form (T-Spline) output
#
# The format is undocumented; it was reverse engineered from one hand-made Quad
# Ball. It is plain text, a half-edge ("link") structure:
#
#   f <a link of this face> 0
#   e <a link of this edge> 1            1 = smooth, 0 = creased
#   v <a link of this vertex> NORTH
#   l <next> <prev> <mate> <vertex> <face> <edge> 0
#   0g <x> <y> <z> <w>                   one per vertex, in v order
#
# Verified against the sample: `next` and `prev` are inverses, `mate` is an
# involution that shares the edge and belongs to the other face, and `vertex` is
# the link's DESTINATION (vert[mate] == vert[prev]). The `next` loop runs
# CLOCKWISE seen from outside on all 24 faces of the sample, which is the
# opposite of the CUBE_FACES convention used everywhere else here, so the loops
# are emitted reversed.
# --------------------------------------------------------------------------- #

TSM_HEADER = """#TS0200

degree 3
cap-type G1CAPS
star-smoothness 0
units 1 meters
end-conditions SUBD_CREASES
star-knot-rule NURCCS"""

TSM_FOOTER = """tol 1.00000000000000008e-05
geom-tol 1.00000000000000008e-05
ver 1044
behavior-version 14.2.0
compat-version 13.3.0"""


def cage_quads(grids):
    """Every cage quad as a 4-tuple of vertex indices, wound +u then +v."""
    out = []
    for g in grids:
        n = g.shape[0] - 1
        for i in range(n):
            for j in range(n):
                out.append((int(g[i, j]), int(g[i + 1, j]),
                            int(g[i + 1, j + 1]), int(g[i, j + 1])))
    return out


def write_tsm(points, faces, scale=1.0):
    """Serialise a polygon cage as a T-Spline control mesh.

    Faces may have ANY arity, not just four. A half-edge structure does not care,
    and the arity restriction was costing real quality: a quad mesher leaves a
    handful of triangles behind, and the only way to be rid of them was to
    midpoint-subdivide the whole mesh, which quadruples the element count and
    lets every bad element breed four more. Letting two triangles through
    instead keeps the good cage.
    """
    # reversed: our faces are counter-clockwise from outside, TSM wants clockwise
    faces = [tuple(f)[::-1] for f in faces]

    edge_id = {}
    for q in faces:
        n = len(q)
        for k in range(n):
            key = (q[k], q[(k + 1) % n])
            key = key if key[0] < key[1] else (key[1], key[0])
            if key not in edge_id:
                edge_id[key] = len(edge_id)

    start = []
    n_links = 0
    for q in faces:
        start.append(n_links)
        n_links += len(q)
    nxt = [0] * n_links
    prv = [0] * n_links
    mate = [-1] * n_links
    dest = [0] * n_links
    fidx = [0] * n_links
    eidx = [0] * n_links
    by_edge = {}
    for f, q in enumerate(faces):
        n = len(q)
        s = start[f]
        for k in range(n):
            i = s + k
            nxt[i] = s + (k + 1) % n
            prv[i] = s + (k - 1) % n
            dest[i] = q[(k + 1) % n]          # link runs q[k] -> q[k+1]
            fidx[i] = f
            key = (q[k], q[(k + 1) % n])
            key = key if key[0] < key[1] else (key[1], key[0])
            eidx[i] = edge_id[key]
            by_edge.setdefault(edge_id[key], []).append(i)
    for e, ls in by_edge.items():
        if len(ls) != 2:
            raise ValueError('edge %d has %d links; the cage is not a closed '
                             'manifold quad mesh' % (e, len(ls)))
        mate[ls[0]], mate[ls[1]] = ls[1], ls[0]

    face_link = list(start)
    edge_link = [by_edge[e][0] for e in range(len(edge_id))]
    vert_link = {}
    for i in range(n_links):
        vert_link.setdefault(dest[i], i)
    n_verts = len(points)
    missing = [v for v in range(n_verts) if v not in vert_link]
    if missing:
        raise ValueError('%d cage vertices are not referenced by any quad'
                         % len(missing))

    out = [TSM_HEADER]
    for i in face_link:
        out.append('f %d 0' % i)
    for i in edge_link:
        out.append('e %d 1' % i)
    for v in range(n_verts):
        out.append('v %d NORTH' % vert_link[v])
    for i in range(n_links):
        out.append('l %d %d %d %d %d %d 0'
                   % (nxt[i], prv[i], mate[i], dest[i], fidx[i], eidx[i]))
    out.append('')
    for v in range(n_verts):
        p = points[v]
        out.append('0g %.17g %.17g %.17g 1'
                   % (p[0] * scale, p[1] * scale, p[2] * scale))
    out.append('')
    out.append(TSM_FOOTER)
    return '\n'.join(out) + '\n'


def limit_positions(points, quads_arr, table):
    """Catmull-Clark limit points of a quad control cage.

    P_limit = (n^2 P + 4 * sum(adjacent) + sum(diagonal)) / (n (n + 5)).
    For a regular vertex that is the familiar (16 P + 4 sum_edge + sum_diag)/36.
    This is what makes it possible to fit the cage so the SURFACE lands on the
    mesh: a subdivision surface does not pass through its own control points.
    """
    a, b, c, d = quads_arr
    valid = table >= 0
    safe = np.where(valid, table, 0)
    adj = (points[safe] * valid[:, :, None]).sum(axis=1)
    n = valid.sum(axis=1).astype(float)

    diag = np.zeros_like(points)
    for src, opp in ((a, c), (b, d), (c, a), (d, b)):
        np.add.at(diag, src, points[opp])

    n2 = n[:, None] ** 2
    return (n2 * points + 4.0 * adj + diag) / (n[:, None] * (n[:, None] + 5.0))


def neighbours_from_quads(n_verts, quads):
    """Ragged neighbour table for a mesh of arbitrary-arity faces."""
    nb = [set() for _ in range(n_verts)]
    for q in quads:
        n = len(q)
        for k in range(n):
            u, v = q[k], q[(k + 1) % n]
            nb[u].add(v)
            nb[v].add(u)
    width = max(len(x) for x in nb)
    tab = np.full((n_verts, width), -1, dtype=np.int64)
    for i, x in enumerate(nb):
        for k, v in enumerate(sorted(x)):
            tab[i, k] = v
    return tab


def subdivision_matrix(n_verts, faces):
    """One Catmull-Clark step as a sparse operator, for faces of ANY arity.

    Returned as a matrix rather than as moved points because the same step has
    to be applied to the cage (to measure) and composed into the limit operator
    (to fit). Doing it once, linearly, keeps those two definitions identical.

    Arity matters because a gmsh cage is not all quads: recombination leaves a
    few hundred triangles on anatomy, and the quad-only version silently dropped
    them. That left holes in the surface being measured and vertices with no
    neighbours at all, which is where the limit formula divided by zero.

    Output faces are always quads - one CC step quadrangulates anything - so
    everything downstream of the first level may assume quads again.

    Returns (S, new_faces), S of shape (n_new, n_verts); new points are
    S @ points. New vertices are ordered originals, then face points in face
    order, then edge points in first-seen edge order.
    """
    from scipy.sparse import coo_matrix
    rows, cols, data = [], [], []

    for f, q in enumerate(faces):
        w = 1.0 / len(q)
        for i in q:
            rows.append(n_verts + f)
            cols.append(i)
            data.append(w)

    edge_faces = {}
    for f, q in enumerate(faces):
        n = len(q)
        for k in range(n):
            u, v = q[k], q[(k + 1) % n]
            key = (u, v) if u < v else (v, u)
            edge_faces.setdefault(key, []).append(f)

    base_e = n_verts + len(faces)
    edge_pt = {}
    for e, (key, fs) in enumerate(edge_faces.items()):
        u, v = key
        r = base_e + e
        edge_pt[key] = r
        w = 1.0 / (2.0 + len(fs))
        rows.extend((r, r))
        cols.extend((u, v))
        data.extend((w, w))
        for f in fs:
            fw = w / len(faces[f])          # the face point, in terms of originals
            for i in faces[f]:
                rows.append(r)
                cols.append(i)
                data.append(fw)

    # the originals move: (F + 2R + (n-3) P) / n, F and R averaged over the
    # incident faces and edges. On a closed manifold #faces == #edges == n.
    val = np.zeros(n_verts)
    for key in edge_faces:
        val[key[0]] += 1
        val[key[1]] += 1
    nv = np.maximum(val, 1.0)

    for q in faces:
        fw = 1.0 / len(q)
        for i in q:
            c = fw / (nv[i] * nv[i])             # F/n, then the outer /n
            for j in q:
                rows.append(i)
                cols.append(j)
                data.append(c)
    for u, v in edge_faces:
        for i in (u, v):
            # 2 * (sum of midpoints / n), then the outer / n; the midpoint
            # itself is (P_u + P_v)/2, so the 2 and the 1/2 cancel
            c = 1.0 / (nv[i] * nv[i])
            rows.extend((i, i))
            cols.extend((u, v))
            data.extend((c, c))
    for i in range(n_verts):
        rows.append(i)
        cols.append(i)
        data.append((nv[i] - 3.0) / nv[i])

    S = coo_matrix((data, (rows, cols)),
                   shape=(base_e + len(edge_faces), n_verts)).tocsr()

    new_faces = []
    for f, q in enumerate(faces):
        fp = n_verts + f
        n = len(q)
        for k in range(n):
            u = q[k]
            pu, nu = q[k - 1], q[(k + 1) % n]
            prev_key = (pu, u) if pu < u else (u, pu)
            next_key = (u, nu) if u < nu else (nu, u)
            new_faces.append((u, edge_pt[next_key], fp, edge_pt[prev_key]))
    return S, new_faces


def catmull_clark(points, faces):
    """One Catmull-Clark subdivision step on a closed mesh of any face arity.

    Needed to measure the Form output honestly: a T-Spline body is the limit
    surface of its control cage, so the cage itself is not what has to match the
    scan. Subdividing twice and taking limit points gives a dense, faithful
    sample of the surface Fusion will actually build.
    """
    S, new_faces = subdivision_matrix(len(points), faces)
    return S @ np.asarray(points, float), new_faces


def limit_surface_mesh(points, quads, levels=2):
    """Dense triangle mesh of the T-Spline limit surface, for measurement."""
    p, q = np.asarray(points, float), [tuple(f) for f in quads]
    if any(len(f) != 4 for f in q):
        levels = max(1, levels)     # one step quadrangulates; the rest assume it
    for _ in range(levels):
        p, q = catmull_clark(p, q)
    tab = neighbours_from_quads(len(p), q)
    qa = (np.array([x[0] for x in q]), np.array([x[1] for x in q]),
          np.array([x[2] for x in q]), np.array([x[3] for x in q]))
    lim = limit_positions(p, qa, tab)
    tris = []
    for a, b, c, d in q:
        tris.append((a, b, c))
        tris.append((a, c, d))
    return lim, np.array(tris, dtype=np.int64)


def limit_matrix(n_verts, quads_arr, table):
    """The Catmull-Clark limit operator as a sparse matrix, so the cage can be
    SOLVED for rather than relaxed toward."""
    from scipy.sparse import coo_matrix
    a, b, c, d = quads_arr
    val = (table >= 0)
    n = val.sum(axis=1).astype(float)
    denom = n * (n + 5.0)

    rows = [np.arange(n_verts)]
    cols = [np.arange(n_verts)]
    data = [n * n / denom]

    for k in range(table.shape[1]):
        m = val[:, k]
        idx = np.flatnonzero(m)
        rows.append(idx)
        cols.append(table[idx, k])
        data.append(4.0 / denom[idx])

    for src, opp in ((a, c), (b, d), (c, a), (d, b)):
        rows.append(src)
        cols.append(opp)
        data.append(1.0 / denom[src])

    M = coo_matrix((np.concatenate(data),
                    (np.concatenate(rows), np.concatenate(cols))),
                   shape=(n_verts, n_verts)).tocsr()
    # Row-normalise. On an all-quad cage every row already sums to 1 and this is
    # a no-op. A cage carrying a few triangles has no opposite corner for the
    # face term, so those rows come out short; renormalising keeps the operator
    # affine -- constants map to constants -- which is what stops the fit
    # translating or shrinking the whole cage.
    from scipy.sparse import diags
    rs = np.asarray(M.sum(axis=1)).ravel()
    rs[rs == 0.0] = 1.0
    return (diags(1.0 / rs) @ M).tocsc()


def limit_operator(n_verts, faces):
    """Catmull-Clark limit operator for a cage of ANY face arity.

    The closed-form mask (n^2 P + 4 sum_adj + sum_diag) / (n(n+5)) holds only on
    an all-quad cage: a triangle has no opposite corner, so the diagonal term is
    simply missing and row-normalising it back to 1 is a guess, not the limit.
    One subdivision step quadrangulates anything, and subdividing does not move
    the limit surface, so the honest operator is the quad mask on the once-refined
    cage composed with the subdivision itself.

    All-quad cages take the direct path unchanged, so the cube-sphere numbers
    stay exactly what they were.
    """
    faces = [tuple(f) for f in faces]
    if all(len(f) == 4 for f in faces):
        table = neighbours_from_quads(n_verts, faces)
        qa = tuple(np.array([f[k] for f in faces], np.int64) for k in range(4))
        return limit_matrix(n_verts, qa, table)

    S, ref = subdivision_matrix(n_verts, faces)
    tab2 = neighbours_from_quads(S.shape[0], ref)
    qa2 = tuple(np.array([f[k] for f in ref], np.int64) for k in range(4))
    Lq = limit_matrix(S.shape[0], qa2, tab2)
    return (Lq[:n_verts, :] @ S).tocsc()


def graph_laplacian(n_verts, table):
    """Uniform graph Laplacian of the cage, used to keep the control net tame."""
    from scipy.sparse import coo_matrix
    val = table >= 0
    n = val.sum(axis=1).astype(float)
    rows = [np.arange(n_verts)]
    cols = [np.arange(n_verts)]
    data = [np.ones(n_verts)]
    for k in range(table.shape[1]):
        m = val[:, k]
        idx = np.flatnonzero(m)
        rows.append(idx)
        cols.append(table[idx, k])
        data.append(-1.0 / n[idx])
    return coo_matrix((np.concatenate(data),
                       (np.concatenate(rows), np.concatenate(cols))),
                      shape=(n_verts, n_verts)).tocsc()


def fit_limit_to_mesh(points, quads_arr, table, proj, iters, lam=1.0e-3, log=None,
                      faces=None):
    """Solve for the control cage whose LIMIT surface sits on the mesh.

    A T-Spline body is the limit surface of its cage, and a subdivision surface
    does not pass through its own control points - a cage lying exactly on the
    scan gives a body several millimetres inside it. So the cage must be solved
    for, not relaxed toward: relaxing diverges here because the limit mask puts
    only n/(n+5) on the diagonal and the rest on the neighbours, leaving Jacobi
    nothing to lean on (the sphere test came out 536% over volume).

    The solve is regularised. Solving L C = T exactly makes the control net
    overshoot wherever the surface is concave, and that multiplied the number of
    buckled quads by five to six - fatal, because Fusion refuses to build a body
    from a self-intersecting T-Spline. Minimising ||L C - T||^2 + lam ||D C||^2
    keeps the net tame; the caller raises lam until no quad is folded.
    """
    from scipy.sparse.linalg import spsolve
    L = (limit_matrix(len(points), quads_arr, table) if faces is None
         else limit_operator(len(points), faces))
    D = graph_laplacian(len(points), table)
    A = (L.T @ L + lam * (D.T @ D)).tocsc()
    for it in range(iters):
        lim = L @ points
        tgt, dist, _f = proj(lim, k=16, exact=False)
        rhs = L.T @ tgt
        points = np.column_stack([spsolve(A, rhs[:, k]) for k in range(3)])
        if log and (it % 5 == 0 or it == iters - 1):
            log('  limit fit %2d/%d  lam %.1e  surface-at-cage max %.4f mean %.4f mm'
                % (it + 1, iters, lam, dist.max() * 10.0, dist.mean() * 10.0))
    lim = L @ points
    _t, dist, _f = proj(lim)
    return points, float(dist.max() * 10.0), float(dist.mean() * 10.0)


def untangle_cage(points, table, quads, passes, log=None):
    """Unbuckle a CONTROL cage with purely local smoothing and no projection.

    The wrapped-cage untangler re-projects onto the mesh, which is exactly wrong
    here: a fitted control net does not lie on the mesh, its limit surface does.
    Smoothing only the folded neighbourhoods leaves the rest of the fit untouched,
    so the accuracy cost stays local instead of global - raising the solve's
    regularisation until folds cleared cost a factor of four in RMS.
    """
    a, b, c, d = quads
    weight = 0.5
    first = None
    for it in range(passes):
        bad, _l = folded_quads(points, quads)
        n_bad = int(bad.sum())
        if first is None:
            first = n_bad
        if n_bad == 0:
            break
        touched = np.zeros(len(points), bool)
        for arr in (a, b, c, d):
            touched[arr[bad]] = True
        for _ring in range(2):
            nb = table[touched]
            touched[nb[nb >= 0]] = True
        points = points + (weight * touched[:, None]) * laplacian(points, table)
        weight = min(0.95, weight * 1.1)
    bad, _l = folded_quads(points, quads)
    if log:
        log('  cage untangle: folded quads %s -> %d' % (first, int(bad.sum())))
    return points, int(bad.sum())


def fit_limit_unfolded(points, quads_arr, table, proj, iters, lam0, passes,
                       log=None, faces=None):
    """Fit the limit surface, then clear any buckle it introduced.

    Solving L C = T exactly makes the control net overshoot wherever the surface
    is concave: measured on the calcaneus it turned 1-21 folded quads into 28-130.
    A light regularisation plus a LOCAL untangle of the result keeps the accuracy
    the fit earned; only if folds survive does the regularisation go up.
    """
    lam = lam0
    best = None
    for attempt in range(5):
        pts, mx, mn = fit_limit_to_mesh(points, quads_arr, table, proj, iters,
                                        lam=lam, log=log if attempt == 0 else None,
                                        faces=faces)
        pts, n_bad = untangle_cage(pts, table, quads_arr, passes, log)
        if best is None or n_bad < best[3]:
            best = (pts, mx, mn, n_bad, lam)
        if n_bad == 0:
            return pts, mx, mn, 0, lam
        if log:
            log('  %d folds survive at lam %.1e; raising it' % (n_bad, lam))
        lam *= 6.0
    return best


def orient_faces(faces):
    """Make every face wind the same way round the surface.

    gmsh orients each classified patch independently, so a cage assembled from
    several patches can have neighbours wound opposite ways. Signed volume then
    means nothing -- measured -46% against the source on a cage that was
    geometrically fine -- and the half-edge mates come out inconsistent, which
    is exactly what a T-Spline reads as a self-intersecting surface.

    BFS the face-adjacency graph and flip any neighbour that traverses a shared
    edge in the SAME direction as the face it came from; two correctly oriented
    faces always traverse their shared edge in opposite directions.
    """
    half = {}
    for fi, f in enumerate(faces):
        n = len(f)
        for k in range(n):
            half.setdefault((f[k], f[(k + 1) % n]), []).append(fi)
    adj = {}
    for (a, b), fs in half.items():
        for fi in fs:
            lst = adj.setdefault(fi, [])
            for fj in half.get((b, a), ()):
                lst.append((fj, False))       # opposite traversal: consistent
            for fj in fs:
                if fj != fi:
                    lst.append((fj, True))    # same traversal: must flip
    seen = [False] * len(faces)
    flip = [False] * len(faces)
    for root in range(len(faces)):
        if seen[root]:
            continue
        seen[root] = True
        stack = [root]
        while stack:
            fi = stack.pop()
            for fj, must in adj.get(fi, ()):
                if not seen[fj]:
                    seen[fj] = True
                    flip[fj] = flip[fi] ^ must
                    stack.append(fj)
    return [tuple(f)[::-1] if flip[i] else tuple(f) for i, f in enumerate(faces)]


def collapse_short_edges(pts, faces, floor, log):
    """Enforce a minimum element edge by welding the ends of shorter edges.

    gmsh's MeshSizeMin governs the point INSERTION field, not what recombination
    then does with the triangles, so the floor is routinely violated: asking for
    1.5 mm produced a 0.30 mm shortest edge. Rather than hope, merge the ends of
    any edge below the floor.

    Union-find over the short edges only, so a chain of them collapses to one
    point instead of leaving a shorter edge behind. A face that loses a corner
    becomes a triangle, which write_tsm now accepts, and a face that collapses
    below three distinct corners is dropped. If the result stops being a closed
    manifold the collapse is abandoned -- a floor is not worth a broken cage.
    """
    if floor <= 0.0:
        return pts, faces
    parent = list(range(len(pts)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    short = 0
    for f in faces:
        n = len(f)
        for k in range(n):
            a, b = f[k], f[(k + 1) % n]
            if np.linalg.norm(pts[a] - pts[b]) < floor:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[max(ra, rb)] = min(ra, rb)
                    short += 1
    if not short:
        return pts, faces

    groups = {}
    for i in range(len(pts)):
        groups.setdefault(find(i), []).append(i)
    order = sorted(groups)
    new_id = {root: i for i, root in enumerate(order)}
    new_pts = np.array([pts[groups[r]].mean(axis=0) for r in order])

    new_faces = []
    for f in faces:
        seq = [new_id[find(v)] for v in f]
        ded = [v for i, v in enumerate(seq) if v != seq[i - 1]]     # cyclic dedup
        if len(ded) >= 3 and len(set(ded)) == len(ded):
            new_faces.append(tuple(ded))

    counts = {}
    for f in new_faces:
        n = len(f)
        for k in range(n):
            a, b = f[k], f[(k + 1) % n]
            key = (a, b) if a < b else (b, a)
            counts[key] = counts.get(key, 0) + 1
    if any(c != 2 for c in counts.values()) or len(new_faces) < 4:
        log('  edge floor %.2f mm abandoned: the collapse broke the manifold'
            % floor)
        return pts, faces
    used = {v for f in new_faces for v in f}
    if len(used) != len(new_pts):
        remap = {v: i for i, v in enumerate(sorted(used))}
        new_pts = new_pts[sorted(used)]
        new_faces = [tuple(remap[v] for v in f) for f in new_faces]
    log('  edge floor %.2f mm: %d short edges collapsed, %d -> %d vertices, '
        '%d -> %d faces'
        % (floor, short, len(pts), len(new_pts), len(faces), len(new_faces)))
    return new_pts, new_faces


def quadrangulate(points, faces):
    """Topological Catmull-Clark step: an n-gon becomes n quads, nothing moves.

    The subdivision rules are deliberately NOT applied. A face point is the
    centroid of its corners and an edge point is the midpoint of its ends, so
    every new vertex lies in the plane of a face it came from - the cage is on
    the surface by construction rather than by projection.
    """
    pts = list(np.asarray(points, float))
    face_pt = []
    for f in faces:
        face_pt.append(len(pts))
        pts.append(np.asarray([points[i] for i in f], float).mean(axis=0))

    edge_pt = {}
    for f in faces:
        n = len(f)
        for k in range(n):
            a, b = int(f[k]), int(f[(k + 1) % n])
            key = (a, b) if a < b else (b, a)
            if key not in edge_pt:
                edge_pt[key] = len(pts)
                pts.append(0.5 * (points[a] + points[b]))

    quads = []
    for i, f in enumerate(faces):
        n = len(f)
        fp = face_pt[i]
        for k in range(n):
            u = int(f[k])
            pu, nu = int(f[k - 1]), int(f[(k + 1) % n])
            pk = (pu, u) if pu < u else (u, pu)
            nk = (u, nu) if u < nu else (nu, u)
            quads.append((u, edge_pt[nk], fp, edge_pt[pk]))
    return np.array(pts), quads


def cage_coverage(points, faces, mesh, samples=40000, seed=20260922):
    """How far the scan is from the cage - the check gmsh's cage needed.

    The fit that follows moves control points; it cannot invent topology. So a
    cage that does not lie over a part of the scan produces a body missing that
    part, and no amount of fitting recovers it. Measuring mesh->cage (not the
    reverse) is what catches that: cage->mesh stays small even when the cage has
    shrunk past a protrusion, because every cage point is still near SOME bit of
    bone.
    """
    tri = []
    for f in faces:
        tri.append((f[0], f[1], f[2]))
        for k in range(3, len(f)):
            tri.append((f[0], f[k - 1], f[k]))
    cage = trimesh.Trimesh(vertices=np.asarray(points, float),
                           faces=np.array(tri, np.int64), process=False)
    # seeded: Auto probes the same cages twice in a session and a wandering
    # answer would read as flakiness rather than as sampling noise
    probe, _fi = trimesh.sample.sample_surface(mesh, int(samples), seed=seed)
    probe = np.vstack([probe, mesh.vertices])
    _c, d, _f = Projector(cage)(probe, k=16, exact=False)
    return {'max_mm': float(d.max() * 10.0),
            'p99_mm': float(np.percentile(d, 99) * 10.0),
            'volume_pct': float(100.0 * (cage.volume - mesh.volume) / mesh.volume)}


def dissolve_doublets(faces):
    """Remove valence-2 vertices: two quads (v,a,x,b) and (v,b,y,a) sharing the
    consecutive edges v-a and v-b become the one quad (a,x,b,y).

    Instant Meshes leaves a few. A doublet is a degenerate T-spline vertex, and
    its two quads share both a and b, so their fan triangulations put the same
    diagonal in four triangles."""
    faces = [list(q) for q in faces]
    removed = 0
    changed = True
    while changed:
        changed = False
        inc = {}
        for fi, q in enumerate(faces):
            if q is None:
                continue
            for v in q:
                inc.setdefault(v, []).append(fi)
        for v, fs in inc.items():
            if len(fs) != 2:
                continue
            f1, f2 = faces[fs[0]], faces[fs[1]]
            if f1 is None or f2 is None or len(f1) != 4 or len(f2) != 4:
                continue
            k = f1.index(v)
            a, x, b = f1[(k + 1) % 4], f1[(k + 2) % 4], f1[(k + 3) % 4]
            j = f2.index(v)
            nb = f2[(j + 1) % 4], f2[(j + 3) % 4]
            if set(nb) != {a, b}:
                continue
            y = f2[(j + 2) % 4]
            if y in (a, b, x):
                continue
            faces[fs[0]] = [a, x, b, y]
            faces[fs[1]] = None
            removed += 1
            changed = True
            break
    return [tuple(q) for q in faces if q is not None], removed


def instant_quad_cage(mesh, job, log):
    """A field-aligned all-quad cage from Instant Meshes (Jakob et al. 2015).
    Points in CENTIMETRES.

    The scan cage quadrangulates triangles, which leaves a valence-3 star at
    every centroid - half the vertices extraordinary. Instant Meshes solves an
    orientation field and a position field on the scan and extracts quads that
    follow them, so the cage is rows of quads with about 2 % extraordinary
    vertices: a clean Form that is easy to edit, and it fits 3x faster. It is
    NOT the more accurate cage. Measured on a 499 mm artificial femur at about
    6000 quads each: scan cage p99 0.12 / max 1.2 mm, this one p99 0.19 / max
    1.8 mm, both with no folds - the scan cage's vertices are scan vertices,
    these are not. Needs `pip install PyNanoInstantMeshes`.

    Its extraction ends in one regular subdivision into pure quads, so the face
    count is about four times the vertex count asked for; the request is sized
    for the final quads to land near hmax_mm.
    """
    import pynanoinstantmeshes as pim

    h = float(job.get('hmax_mm', 8.0))
    cap = int(job.get('cage_quad_cap', 24000))
    area_mm2 = float(mesh.area) * 100.0
    quads = max(200, min(int(area_mm2 / (h * h)), cap))
    v = (np.asarray(mesh.vertices, float) * 10.0).astype(np.float32)
    f = np.asarray(mesh.faces).astype(np.uint32)
    V, F = pim.remesh(v, f, max(50, quads // 4), rosy=4, posy=4,
                      smooth_iter=int(job.get('instant_smooth', 2)),
                      deterministic=True)
    V = np.asarray(V, float)
    faces = []
    for q in np.asarray(F):
        u = []
        for i in (int(x) for x in q):
            if not u or u[-1] != i:          # a triangle comes back with a repeated corner
                u.append(i)
        if len(u) > 1 and u[0] == u[-1]:
            u.pop()
        if len(u) >= 3:
            faces.append(tuple(u))
    faces, n_doublets = dissolve_doublets(faces)
    used = np.unique(np.concatenate([np.array(q) for q in faces]))
    remap = -np.ones(len(V), np.int64)
    remap[used] = np.arange(len(used))
    V = V[used]
    faces = [tuple(int(remap[i]) for i in q) for q in faces]

    # closed means every POLYGON edge is used twice; a fan triangulation is no
    # test, because its diagonals are not cage edges
    uses = {}
    for q in faces:
        for k in range(len(q)):
            e = (min(q[k], q[k - 1]), max(q[k], q[k - 1]))
            uses[e] = uses.get(e, 0) + 1
    if any(c != 2 for c in uses.values()):
        raise RuntimeError('Instant Meshes returned an open cage')
    faces = orient_faces(faces)
    tri = [(q[0], q[k - 1], q[k]) for q in faces for k in range(2, len(q))]
    if trimesh.Trimesh(vertices=V, faces=np.array(tri), process=False).volume < 0:
        faces = [tuple(reversed(q)) for q in faces]
    if n_doublets:
        log('  instant cage: %d doublet vertices dissolved' % n_doublets)
    val = np.bincount(np.concatenate([np.array(q) for q in faces]), minlength=len(V))
    nq = sum(1 for q in faces if len(q) == 4)
    log('  instant cage: %d quads + %d other faces, %d vertices, %.1f %% '
        'extraordinary (element %.2f mm asked)'
        % (nq, len(faces) - nq, len(V), 100.0 * float((val != 4).mean()), h))
    return V * 0.1, faces


def make_cage(mesh, job, log):
    """The cage for job['cage_source']; Instant Meshes falls back to the scan cage."""
    if job.get('cage_source') == 'instant':
        try:
            return instant_quad_cage(mesh, job, log), 'instant'
        except Exception as exc:
            log('  instant cage unavailable (%s: %s), using the scan cage'
                % (type(exc).__name__, exc))
    return scan_quad_cage(mesh, job, log), 'scan'


def scan_quad_cage(mesh, job, log):
    """An all-quad cage taken straight off the scan. Points in CENTIMETRES.

    gmsh's cage comes from a global reparametrisation of the scan, and on a bone
    that parametrisation is so distorted that the cage wanders off the surface -
    measured here at 6 mm away and 11.8% short on volume, with the scan up to
    11.3 mm from the cage it was supposed to sit on. Asking for finer elements
    made it worse, not better, because a finer mesh in a distorted parameter
    space just resolves the distortion.

    This route never leaves the scan: decimate it to the triangle count that
    gives the requested element size, then quadrangulate. Each triangle becomes
    three quads whose vertices are its corners, edge midpoints and centroid, all
    on the decimated surface. Density is a COUNT, so it is exact and monotone,
    and the cage is watertight whenever the decimation is.

    The trade is element shape: a quadrangulated triangle leaves a valence-3 star
    at every centroid, where gmsh's recombination gives cleaner rows. T-Splines
    take star points in their stride, and the cage may be relaxed afterwards.
    """
    import pyvista as pv

    h = float(job.get('hmax_mm', 8.0))
    cap = int(job.get('cage_quad_cap', 24000))
    area_mm2 = float(mesh.area) * 100.0

    # a triangle of edge L quadrangulates into quads of edge about L/2, and
    # F triangles over area A have L = sqrt(4A / (sqrt(3) F)); solving for the
    # count that lands the QUAD edge on h. Measured against the achieved median
    # edge on this scan it is within 6%.
    target = int(round(area_mm2 / (math.sqrt(3.0) * h * h)))
    target = max(40, min(target, cap // 3))

    v_mm = np.asarray(mesh.vertices, float) * 10.0
    f = np.asarray(mesh.faces)
    pd = pv.PolyData(v_mm, np.hstack([np.full((len(f), 1), 3), f]).ravel())
    if target < len(f):
        dec = pd.decimate(1.0 - float(target) / len(f), volume_preservation=True)
        dec = dec.clean().triangulate()
    else:
        dec = pd.clean().triangulate()
    vd = np.asarray(dec.points, float)
    fd = np.asarray(dec.faces).reshape(-1, 4)[:, 1:]

    base = trimesh.Trimesh(vertices=vd, faces=fd, process=True)
    if not base.is_watertight:
        base.merge_vertices()
    if not base.is_watertight:
        raise RuntimeError('decimating the scan to %d triangles left it open; '
                           'ask for a larger element' % target)
    base.fix_normals()
    vd = np.asarray(base.vertices, float)
    fd = [tuple(int(x) for x in t) for t in np.asarray(base.faces)]

    pts, quads = quadrangulate(vd, fd)

    # A quadrangulated triangle inherits its shape: in an obtuse triangle the
    # centroid sits close to the long edge's midpoint, and the quad between them
    # comes out 45:1. The decimated TRIANGLES are fine (ratio max 7.7, shortest
    # edge 0.9 mm) - it is the quadrangulation that makes the slivers, so this is
    # the place to collapse them. A quad that loses a corner becomes a triangle,
    # which the TSM writer and the limit operator both take.
    floor = float(job.get('cage_edge_floor', 0.25)) * h
    if floor > 0.0:
        pts, quads = collapse_short_edges(pts, quads, floor, log)

    quads = orient_faces(quads)
    tri = []
    for q in quads:
        for k in range(2, len(q)):
            tri.append((q[0], q[k - 1], q[k]))
    if trimesh.Trimesh(vertices=pts, faces=np.array(tri),
                       process=False).volume < 0:
        quads = [tuple(reversed(q)) for q in quads]

    qs = np.array([q for q in quads if len(q) == 4], np.int64)
    L = edge_lengths(pts, (qs[:, 0], qs[:, 1], qs[:, 2], qs[:, 3]))
    nq = len(qs)
    log('  scan cage: %d quads + %d triangles from %d decimated triangles, '
        '%d vertices (element %.2f mm asked, %.2f mm median)'
        % (nq, len(quads) - nq, len(fd), len(pts), h, float(np.median(L))))
    return pts * 0.1, quads          # mm -> cm


def gmsh_cage(mesh, job, log):
    """A quad cage from gmsh, already sitting on the mesh. Points in CENTIMETRES.

    This replaces cube_sphere for the quad-cage method. The difference that
    matters is not the mesher but the inversion: element SIZE is the input here
    and the element count follows, where a cube-sphere fixes the count first and
    lets size fall out wherever the surface happens to need area. That is why a
    fixed grid cannot hold an edge-length ratio on anatomy and this can.
    """
    import os
    import tempfile
    import gmsh
    import pyvista as pv

    hmax = float(job.get('hmax_mm', 8.0))
    hmin = float(job.get('hmin_mm', 1.5))
    curv = float(job.get('curvature', 16.0))
    budget = int(job.get('cage_input_tris', 20000))

    v_mm = np.asarray(mesh.vertices, float) * 10.0
    f = np.asarray(mesh.faces)
    if len(f) > budget:
        pd = pv.PolyData(v_mm, np.hstack([np.full((len(f), 1), 3), f]).ravel())
        sm = pd.decimate(1.0 - float(budget) / len(f)).clean().triangulate()
        v_in = np.asarray(sm.points, float)
        f_in = np.asarray(sm.faces).reshape(-1, 4)[:, 1:]
        log('  cage input decimated %d -> %d triangles' % (len(f), len(f_in)))
    else:
        v_in, f_in = v_mm, f

    tmp = os.path.join(tempfile.gettempdir(), 'mtp_cage_%d.stl' % os.getpid())
    trimesh.Trimesh(vertices=v_in, faces=f_in, process=True).export(tmp)
    try:
        gmsh.initialize()
        gmsh.option.setNumber('General.Terminal', 0)
        gmsh.merge(tmp)
        # pi, not a feature angle: on a noisy scan a 60 degree angle shatters
        # the surface into thousands of patches and classifySurfaces never
        # returns. At pi it is instant and the reparametrisation still splits
        # wherever it has to.
        gmsh.model.mesh.classifySurfaces(math.pi, True, True, math.pi)
        gmsh.model.mesh.createGeometry()
        gmsh.model.geo.synchronize()
        for k, val in (('Mesh.MeshSizeMin', hmin), ('Mesh.MeshSizeMax', hmax),
                       ('Mesh.MeshSizeFromCurvature', curv),
                       ('Mesh.MeshSizeExtendFromBoundary', 0),
                       ('Mesh.Algorithm', 8), ('Mesh.RecombineAll', 1),
                       ('Mesh.RecombinationAlgorithm', 3),
                       ('Mesh.SubdivisionAlgorithm', 0)):
            gmsh.option.setNumber(k, val)
        gmsh.model.mesh.generate(2)
        tags, coords, _p = gmsh.model.mesh.getNodes()
        pts = np.asarray(coords, float).reshape(-1, 3)
        remap = {int(t): i for i, t in enumerate(tags)}
        etypes, _eids, enodes = gmsh.model.mesh.getElements(2)
        faces = []
        for t, nodes in zip(etypes, enodes):
            arity = {2: 3, 3: 4}.get(int(t))
            if not arity:
                continue
            arr = np.array([remap[int(x)] for x in nodes], np.int64).reshape(-1, arity)
            faces.extend(tuple(int(x) for x in row) for row in arr)
    finally:
        try:
            gmsh.finalize()
        except Exception:
            pass
        try:
            os.remove(tmp)
        except Exception:
            pass

    if not faces:
        raise RuntimeError('gmsh produced no surface elements for the cage')

    used = sorted({v for fc in faces for v in fc})
    idx = {v: i for i, v in enumerate(used)}
    pts = pts[used]
    faces = [tuple(idx[v] for v in fc) for fc in faces]
    faces = orient_faces(faces)
    pts, faces = collapse_short_edges(pts, faces, hmin * float(job.get('edge_floor', 0.6)), log)

    tri = []
    for fc in faces:
        tri.append((fc[0], fc[1], fc[2]))
        if len(fc) == 4:
            tri.append((fc[0], fc[2], fc[3]))
    if trimesh.Trimesh(vertices=pts, faces=np.array(tri), process=False).volume < 0:
        faces = [tuple(reversed(fc)) for fc in faces]

    nq = sum(1 for fc in faces if len(fc) == 4)
    log('  gmsh cage: %d quads + %d triangles, %d vertices (size %.1f-%.1f mm, '
        'curvature %d)' % (nq, len(faces) - nq, len(pts), hmin, hmax, int(curv)))
    return pts * 0.1, faces          # mm -> cm


def probe_cages(job, log):
    """Build cages at several element sizes and measure them, without fitting.

    The fit is the expensive part - 8 to 27 s on a distal radius - and the cage
    is not: 0.2 s once the imports are warm. So the way to choose an element size
    is to build the cages and look at them, rather than to guess from a formula
    and then discover the answer twenty seconds later.

    What is reported is mesh->cage, because that is the half that says whether
    the cage has lost a feature; see cage_coverage. Measured against the fitted
    result on a distal radius over 8.0 to 1.2 mm, the final worst-of-both-
    directions deviation came out at 0.44 to 0.53 times the cage's own max for
    every size where the fit behaved, so half the cage max is the working
    prediction and the cage max itself is the conservative one.
    """
    coords = np.array(job['coords'], float).reshape(-1, 3)
    faces = np.array(job['idx'], np.int64).reshape(-1, 3)
    mesh = trimesh.Trimesh(vertices=coords, faces=faces, process=True)
    diag = float(np.linalg.norm(mesh.bounds[1] - mesh.bounds[0])) * 10.0
    out = []
    for h in job['probe']:
        row = {'element_mm': float(h)}
        try:
            (pts, cf), used = make_cage(mesh, dict(job, hmax_mm=float(h)),
                                        lambda m: None)
            cov = cage_coverage(pts, cf, mesh,
                                samples=int(job.get('probe_samples', 20000)))
            row.update(cov)
            row['quads'] = int(sum(1 for f in cf if len(f) == 4))
            row['vertices'] = int(len(pts))
            # what the fit is expected to land on, and the pessimistic reading
            # the fit lands at 0.44-0.53 of a scan cage's max (distal radius) but
            # only 0.76-0.86 of an instant cage's (proxy femur): the instant
            # cage is not a subset of the scan's vertices, so it starts further off
            row['cage_source'] = used
            row['predicted_mm'] = (0.9 if used == 'instant' else 0.5) * cov['max_mm']
            row['pessimistic_mm'] = cov['max_mm']
        except Exception as exc:
            row['error'] = '%s: %s' % (type(exc).__name__, exc)
        out.append(row)
        log('  element %5.2f mm -> %s' % (h, row.get('error') or
            'cage max %.3f p99 %.3f mm, %d quads'
            % (row['max_mm'], row['p99_mm'], row['quads'])))
    return {'probe': out, 'diagonal_mm': diag,
            'triangles': int(len(mesh.faces)),
            'watertight': bool(mesh.is_watertight),
            'area_mm2': float(mesh.area) * 100.0,
            'volume_mm3': float(mesh.volume) * 1000.0}


def build(job, log):
    t0 = time.time()
    coords = np.array(job['coords'], float).reshape(-1, 3)
    faces = np.array(job['idx'], np.int64).reshape(-1, 3)
    mesh = trimesh.Trimesh(vertices=coords, faces=faces, process=True)
    log('mesh %d verts %d faces watertight=%s'
        % (len(mesh.vertices), len(mesh.faces), mesh.is_watertight))

    n = int(job.get('cage', 24))
    m = n + 1
    cage_faces = None
    grids = None
    cage_source_used = None
    cage_cover = None

    if job.get('cage_source') in ('gmsh', 'scan', 'instant'):
        if job.get('output') != 'tsm':
            raise RuntimeError('the quad cage is a Form (T-Spline) cage; '
                               'the NURBS path needs the cube-sphere')
        source = job['cage_source']
        h = float(job.get('hmax_mm', 8.0))
        if source == 'gmsh':
            points, cage_faces = gmsh_cage(mesh, job, log)
            cover = cage_coverage(points, cage_faces, mesh)
            log('  cage covers the scan to max %.3f p99 %.3f mm, volume %+.2f %%'
                % (cover['max_mm'], cover['p99_mm'], cover['volume_pct']))
            # a control cage that does not lie OVER a part of the scan cannot be
            # fitted onto it later - the fit moves points, not topology. One
            # element of slack is generous; beyond that the cage has lost
            # something and a finer one loses more, so change route instead.
            if (cover['p99_mm'] > float(job.get('coverage_slack', 1.0)) * h
                    and job.get('cage_fallback', True)):
                log('  that is more than one element away: gmsh\'s '
                    'reparametrisation has lost part of the scan, '
                    'falling back to the scan-derived cage')
                source = 'scan'
        if source in ('scan', 'instant'):
            (points, cage_faces), source = make_cage(mesh, dict(job, cage_source=source), log)
            cover = cage_coverage(points, cage_faces, mesh)
            log('  cage covers the scan to max %.3f p99 %.3f mm, volume %+.2f %%'
                % (cover['max_mm'], cover['p99_mm'], cover['volume_pct']))
        cage_source_used = source
        cage_cover = cover
        proj = Projector(mesh)
        table = neighbours_from_quads(len(points), cage_faces)
        qs = np.array([f for f in cage_faces if len(f) == 4], np.int64)
        qa = (qs[:, 0].copy(), qs[:, 1].copy(), qs[:, 2].copy(), qs[:, 3].copy())
        quads = [list(f) for f in cage_faces if len(f) == 4]
        cage_folds = 0
        limit = float(job.get('aspect_limit', 0.0)) or 2.0
        cage_quality = {'before': aspect_report(points, qa, limit)}
        # No shrink wrap: the cage is already on the scan, which is the point of
        # generating it there instead of deforming a sphere onto it. A shape
        # relaxation is still worth having - a quadrangulated triangle leaves
        # elongated elements - but it is tangential and fold-safe, so the cage
        # stays on the surface.
        n_sq = int(job.get('cage_square_iters', 0))
        if n_sq > 0:
            _c, _d, fid = proj(points, k=16, exact=False)
            points, fid = relax_square(
                points, table, qa, proj, fid, n_sq,
                float(job.get('square_weight', 0.6)),
                float(job.get('min_edge_mm', 0.0)) * 0.1, limit, log,
                step_cap=float(job.get('square_step_cap', 0.10)))
        cage_quality['after'] = aspect_report(points, qa, limit)
    else:
        verts, grids = cube_sphere(n)
        table = cage_neighbours(len(verts), grids)

        centre = mesh.bounding_box.centroid
        radius = float(np.linalg.norm(mesh.vertices - centre, axis=1).max()) * 1.05
        points = verts * radius + centre
        log('cage %d verts, start sphere r=%.2f mm' % (len(points), radius * 10.0))

        proj = Projector(mesh)
        points, cage_folds, cage_quality = shrink_wrap(
            mesh, proj, points, table, grids,
            iters=int(job.get('wrap_iters', 50)),
            smooth_start=float(job.get('smooth_start', 0.55)),
            smooth_end=float(job.get('smooth_end', 0.06)),
            sample_count=int(job.get('samples', 120000)),
            tangential_iters=int(job.get('relax_iters', 60)),
            untangle_passes=int(job.get('untangle', 12)),
            log=log,
            aspect_limit=float(job.get('aspect_limit', 0.0)),
            aspect_iters=int(job.get('aspect_iters', 0)),
            aspect_gain=float(job.get('aspect_gain', 1.2)),
            square_iters=int(job.get('square_iters', 0)),
            square_weight=float(job.get('square_weight', 0.5)),
            min_edge=float(job.get('min_edge_mm', 0.0)) * 0.1,
            square_step_cap=float(job.get('square_step_cap', 0.10)))

    _c, dcage, _f = proj(points)
    log('cage sits on the mesh to %.5f mm' % (dcage.max() * 10.0))

    # ---- Form / T-Spline output -------------------------------------------- #
    if job.get('output') == 'tsm':
        if grids is not None:
            quads = cage_quads(grids)
            qa = quad_arrays(grids)
        pts, lim_max, lim_mean, lim_folds, lam_used = fit_limit_unfolded(
            points, qa, table, proj, int(job.get('limit_iters', 20)),
            float(job.get('limit_fair', 1.0e-3)),
            int(job.get('cage_untangle', 60)), log, faces=cage_faces)
        log('control cage fitted so the LIMIT surface lands on the mesh')

        # the FULL face list, not the quad subset: dropping the triangles
        # left holes in the surface being measured and vertices with no
        # neighbours at all, which is where the limit formula produced NaN
        lim_pts, lim_tris = limit_surface_mesh(
            pts, cage_faces if cage_faces is not None else quads,
            int(job.get('eval_levels', 2)))
        _c, d_out, _f = proj(lim_pts)
        surf = trimesh.Trimesh(vertices=lim_pts, faces=lim_tris, process=False)
        proj_s = Projector(surf)
        probe, _fi = trimesh.sample.sample_surface(mesh, int(job.get('probe', 60000)))
        probe = np.vstack([probe, mesh.vertices])
        _c2, d_back, _f2 = proj_s(probe)
        folded, _l = folded_quads(pts, qa)


        return {
            'tsm': write_tsm(pts, cage_faces if cage_faces is not None else quads),
            'cage_points': pts.tolist(),
            'stats': {
                'output': 'tsm',
                'cage_source': cage_source_used,
                'cage_cover': cage_cover,
                'cage_square_iters': int(job.get('cage_square_iters', 0)),
                'cage_element_mm': (float(job.get('hmax_mm', 8.0))
                                    if cage_source_used else None),
                'cage': n,
                'cage_verts': int(len(pts)),
                'cage_quads': int(len(quads)),
                'cage_folded_quads': int(folded.sum()),
                # cage_quality is measured on the relaxed cage; the Form path
                # then re-solves it so its limit surface lands on the mesh,
                # which moves every point. Measure again on the cage that is
                # actually emitted, or the reported shape is not the one the
                # user will edit in the Form environment.
                'cage_quality': (
                    dict(cage_quality,
                         emitted=aspect_report(pts, qa,
                                               float(job.get('aspect_limit', 3.0))))
                    if cage_quality else
                    {'emitted': aspect_report(pts, qa,
                                              float(job.get('aspect_limit', 3.0)))}),
                'limit_regularisation': lam_used,
                'limit_at_cage_max_mm': lim_max,
                'limit_at_cage_mean_mm': lim_mean,
                'surface_tris': int(len(lim_tris)),
                'probe_points': int(len(probe)),
                'signed_volume_mm3': float(surf.volume) * 1000.0,
                'mesh_volume_mm3': float(mesh.volume) * 1000.0,
                'winding_consistent': bool(surf.is_winding_consistent),
                'surf_to_mesh_max_mm': float(d_out.max() * 10.0),
                'surf_to_mesh_p99_mm': float(np.percentile(d_out, 99) * 10.0),
                'surf_to_mesh_rms_mm': float(np.sqrt((d_out ** 2).mean()) * 10.0),
                'mesh_to_surf_max_mm': float(d_back.max() * 10.0),
                'mesh_to_surf_p99_mm': float(np.percentile(d_back, 99) * 10.0),
                'mesh_to_surf_rms_mm': float(np.sqrt((d_back ** 2).mean()) * 10.0),
                'within_0p5_pct': float((d_out * 10.0 < 0.5).mean() * 100.0),
                'within_1_pct': float((d_out * 10.0 < 1.0).mean() * 100.0),
                'within_2_pct': float((d_out * 10.0 < 2.0).mean() * 100.0),
                'secs': round(time.time() - t0, 2),
            },
        }

    P = int(job.get('ctrl', max(8, int(round(0.7 * m)))))
    P = max(6, min(P, m))
    lam = float(job.get('fair', 2.0e-4))
    B, kv = design(P, m)
    log('fit: %d data / %d control points per side, fairness %.2e' % (m, P, lam))

    # every cube edge is fitted ONCE, from the shared row of cage points, so the
    # two patches that meet there get the identical boundary curve
    edge_cache = {}

    def edge_controls(ids, data):
        a, b = int(ids[0]), int(ids[-1])
        flip = a > b
        key = (b, a) if flip else (a, b)
        if key not in edge_cache:
            edge_cache[key] = fit_curve_fixed_ends(
                data[::-1] if flip else data, B, lam)
        c = edge_cache[key]
        return c[::-1] if flip else c

    patch_nets = []
    face_loops = []
    for f, g in enumerate(grids):
        D = points[g]                                     # (m, m, 3) as [u][v]
        C = fit_patch_constrained(
            D, B, lam,
            edge_controls(g[0, :], D[0, :, :]),
            edge_controls(g[-1, :], D[-1, :, :]),
            edge_controls(g[:, 0], D[:, 0, :]),
            edge_controls(g[:, -1], D[:, -1, :]))
        patch_nets.append(C)

        # Loop counter-clockwise in (u, v). u x v is the outward normal by the way
        # CUBE_FACES is built, so this loop is counter-clockwise seen from outside
        # and no face needs its parameterisation reversed.
        walk = [
            (int(g[0, 0]), int(g[-1, 0])),        # v=0, u rising
            (int(g[-1, 0]), int(g[-1, -1])),      # u=1, v rising
            (int(g[-1, -1]), int(g[0, -1])),      # v=1, u falling
            (int(g[0, -1]), int(g[0, 0])),        # u=0, v falling
        ]
        loop = []
        for a, b in walk:
            flip = a > b
            key = '%d-%d' % ((b, a) if flip else (a, b))
            loop.append({'edge': key, 'reversed': bool(flip)})
        face_loops.append(loop)

    corners = sorted({int(v) for g in grids
                      for v in (g[0, 0], g[-1, 0], g[-1, -1], g[0, -1])})
    topo = {
        'vertices': {str(c): (points[c]).tolist() for c in corners},
        'edges': {('%d-%d' % k): {'controlPoints': v.tolist(),
                                  'knots': kv.tolist(), 'degree': 3,
                                  'start': k[0], 'end': k[1]}
                  for k, v in edge_cache.items()},
        'faces': face_loops,
    }

    # Measure both ways against a TESSELLATION of the fitted surface, not against
    # a point sample of it: nearest-sample distance charges the fit for the
    # sampling gaps, which on a part this size was worth several millimetres.
    grid_res = int(job.get('eval_res', 4 * n))
    verts_s, tris_s = [], []
    base = 0
    for C in patch_nets:
        e = eval_surface(C, kv, kv, grid_res, grid_res)
        verts_s.append(e.reshape(-1, 3))
        r = grid_res
        ii, jj = np.meshgrid(np.arange(r - 1), np.arange(r - 1), indexing='ij')
        a = (base + ii * r + jj).ravel()
        b = (base + ii * r + jj + 1).ravel()
        c = (base + (ii + 1) * r + jj + 1).ravel()
        d = (base + (ii + 1) * r + jj).ravel()
        # a=(u,v)  b=(u,v+1)  c=(u+1,v+1)  d=(u+1,v).  Winding a-d-c / a-c-b gives
        # the normal along +u x +v, which CUBE_FACES defines as outward.
        tris_s.append(np.concatenate([np.stack([a, d, c], 1),
                                      np.stack([a, c, b], 1)]))
        base += r * r
    verts_s = np.vstack(verts_s)
    tris_s = np.vstack(tris_s)

    # Orientation check. Comparing surface facets against the mesh normal was
    # tried and dropped: on a creased scan it flags facets that are perfectly fine,
    # the same false positive the cage fold test had. Signed volume against the
    # mesh volume is the honest global statement, and the intrinsic cage fold count
    # is the honest local one.
    surf_mesh_tmp = trimesh.Trimesh(vertices=verts_s, faces=tris_s, process=False)
    signed_vol = float(surf_mesh_tmp.volume)

    _c, d_out, _f = proj(verts_s)                          # surface -> mesh
    proj_s = Projector(trimesh.Trimesh(vertices=verts_s, faces=tris_s,
                                       process=False))
    probe, _fi = trimesh.sample.sample_surface(mesh, int(job.get('probe', 60000)))
    probe = np.vstack([probe, mesh.vertices])
    _c2, d_back, _f2 = proj_s(probe)                       # mesh -> surface

    patches = []
    for C in patch_nets:
        patches.append({
            'controlPoints': C.reshape(-1, 3).tolist(),    # [u][v], u major
            'countU': int(C.shape[0]),
            'countV': int(C.shape[1]),
            'knotsU': kv.tolist(),
            'knotsV': kv.tolist(),
            'degree': 3,
        })

    result = {
        'patches': patches,
        'topology': topo,
        'grids': [g.tolist() for g in grids],
        'cage_points': points.tolist(),
        'stats': {
            'cage': n,
            'cage_verts': int(len(points)),
            'control_per_side': P,
            'fairness': lam,
            'cage_on_mesh_max_mm': float(dcage.max() * 10.0),
            'cage_folded_quads': int(cage_folds),
            'cage_quality': cage_quality,
            'signed_volume_mm3': signed_vol * 1000.0,
            'mesh_volume_mm3': float(mesh.volume) * 1000.0,
            'winding_consistent': bool(surf_mesh_tmp.is_winding_consistent),
            'surface_tris': int(len(tris_s)),
            'probe_points': int(len(probe)),
            'surf_to_mesh_max_mm': float(d_out.max() * 10.0),
            'surf_to_mesh_p99_mm': float(np.percentile(d_out, 99) * 10.0),
            'surf_to_mesh_rms_mm': float(np.sqrt((d_out ** 2).mean()) * 10.0),
            'mesh_to_surf_max_mm': float(d_back.max() * 10.0),
            'mesh_to_surf_p99_mm': float(np.percentile(d_back, 99) * 10.0),
            'mesh_to_surf_rms_mm': float(np.sqrt((d_back ** 2).mean()) * 10.0),
            'within_0p5_pct': float((d_out * 10.0 < 0.5).mean() * 100.0),
            'within_1_pct': float((d_out * 10.0 < 1.0).mean() * 100.0),
            'within_2_pct': float((d_out * 10.0 < 2.0).mean() * 100.0),
            'worst_at_mm': (np.asarray(verts_s[int(d_out.argmax())]) * 10.0).tolist(),
            'secs': round(time.time() - t0, 2),
        },
    }
    return result


def main():
    job = json.load(open(sys.argv[1], encoding='utf-8'))
    lines = []

    def log(msg):
        lines.append(msg)
        print(msg, flush=True)

    try:
        res = probe_cages(job, log) if job.get('probe') else build(job, log)
        res['log'] = lines
        res['ok'] = True
    except Exception:
        import traceback
        res = {'ok': False, 'error': traceback.format_exc(), 'log': lines}
    with open(sys.argv[2], 'w', encoding='utf-8') as fh:
        json.dump(res, fh)


if __name__ == '__main__':
    main()
