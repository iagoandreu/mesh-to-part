"""Generate the toolbar icon for the Wrap (offset) command, from the real solids.

    py -3.10 tools/make_wrap_icon.py

WHY IT IS RENDERED AND NOT DRAWN
    Every earlier version of this icon was drawn by hand from bezier controls,
    and every one came out as something else: a hat, a ribbon, a sling, a box
    standing on a tray. The shape is a wrap - a body grown by an offset, over the
    body it came from - and the honest way to draw that is to project the two
    solids that actually demonstrate it.

    Their geometry is exported from Fusion once into wrap_icon_geometry.json and
    kept beside this script, so the icon rebuilds without Fusion running.

THE FUSION STYLE, MEASURED OFF A SCREENSHOT OF THE REAL RIBBON
    Fusion's own icons carry NO outline. They are flat-shaded solids, light on
    the dark ribbon, and all the separation comes from shading: a light top face
    against a darker side. Two materials carry the meaning - grey for the
    geometry you started with, blue for what the operation produces.

    Sampled from that screenshot: ribbon #3B4453, greys #D9D9D9 / #BBBBBB, blues
    #B2DFFF / #89CBFA / #67B1E6.

    This keeps that palette and adds a thick outline in the ribbon's own dark
    grey. Without one the icon works only on a dark ribbon; with it, it reads on
    both themes, which an add-in has no way to choose between.
"""
import json
import math
import os

import numpy as np
from PIL import Image, ImageDraw, ImageFilter

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), 'resources', 'Wrap')
GEOM = os.path.join(HERE, 'wrap_icon_geometry.json')

INK = (59, 68, 83, 255)            # #3B4453, the ribbon's own dark grey
BASE = ((217, 217, 217, 255), (187, 187, 187, 255), (154, 158, 165, 255))
WRAP = ((178, 223, 255, 255), (137, 203, 250, 255), (103, 177, 230, 255))

# the light comes from up and to the left, as it does in every Fusion icon
LIGHT = np.array([-0.35, -0.55, 0.76])
LIGHT /= np.linalg.norm(LIGHT)


def load():
    with open(GEOM, encoding='utf-8') as fh:
        d = json.load(fh)
    return [(p['name'], np.array(p['coords'], float).reshape(-1, 3),
             np.array(p['idx'], np.int64).reshape(-1, 3)) for p in d['parts']]


def iso_matrix(turn_deg=45.0):
    """Isometric: turn about z, then tip 35.264 degrees down."""
    a = math.radians(turn_deg)
    b = math.atan(1.0 / math.sqrt(2.0))
    rz = np.array([[math.cos(a), -math.sin(a), 0],
                   [math.sin(a), math.cos(a), 0], [0, 0, 1]], float)
    rx = np.array([[1, 0, 0], [0, math.cos(b), -math.sin(b)],
                   [0, math.sin(b), math.cos(b)]], float)
    return rx @ rz


def tone(normal, tones):
    """Flat shading quantised to three steps, so coplanar triangles match.

    Quantising matters more than it sounds. A smooth ramp makes the tessellation
    of the wrap's rounded edges show up as banding by 64 px, and Fusion's own
    icons use two or three flat tones per material and nothing between.
    """
    d = float(np.dot(normal, LIGHT))
    return tones[0] if d > 0.62 else (tones[1] if d > 0.12 else tones[2])


# Chosen by measuring, not by taste: sweeping the turn in 10 degree steps and
# counting pixels, the share of the original cube still visible runs from 8 % at
# the front to 33 % at the back. Under about 20 % the grey is gone by 16 px and
# the icon is a blue blob; 170 degrees leaves 31 %, so both materials still read
# at the smallest size, and it is still an isometric three-quarter view rather
# than the flat elevation the extreme gives.
TURN = 170.0


def render(size, supersample=8, edge_frac=0.055, feature_deg=25.0,
           turn_deg=TURN):
    S = size * supersample
    R = iso_matrix(turn_deg)
    parts = load()

    allv = np.vstack([v for _n, v, _f in parts])
    # the centre of gravity of the pair goes to the origin before the shot, so
    # the framing does not depend on where in the model the solids happened to be
    centre = 0.5 * (allv.min(axis=0) + allv.max(axis=0))
    proj = [(n, (v - centre) @ R.T, f) for n, v, f in parts]

    pts = np.vstack([p[:, :2] for _n, p, _f in proj])
    span = max(pts[:, 0].max() - pts[:, 0].min(),
               pts[:, 1].max() - pts[:, 1].min())
    margin = 0.075
    scale = S * (1.0 - 2.0 * margin) / span
    mid = 0.5 * (pts.min(axis=0) + pts.max(axis=0))

    def to_px(p):
        return np.column_stack([(p[:, 0] - mid[0]) * scale + S * 0.5,
                                S * 0.5 - (p[:, 1] - mid[1]) * scale])

    ew = max(1, int(round(S * edge_frac)))

    # A Z-BUFFER, not a painter's sort. The Cube is twelve triangles, so each one
    # is half a 30 mm face, and sorting by centroid depth draws whole grey
    # triangles in front of blue ones they are actually behind - it showed up as
    # grey wedges lying across the wrap. Centroid order is only ever a proxy;
    # per-pixel depth is the answer.
    zbuf = np.full((S, S), -1.0e30)
    rgb = np.zeros((S, S, 4), np.uint8)
    creases = []

    for name, p, f in proj:
        tones = WRAP if 'wrap' in name.lower() else BASE
        xy = to_px(p)
        a, b, c = p[f[:, 0]], p[f[:, 1]], p[f[:, 2]]
        n = np.cross(b - a, c - a)
        ln = np.linalg.norm(n, axis=1)
        n = n / np.maximum(ln, 1e-30)[:, None]
        z = p[:, 2]

        for i in np.nonzero((ln > 1e-12) & (n[:, 2] > 0.0))[0]:
            i0, i1, i2 = f[i]
            x0, y0 = xy[i0]
            x1, y1 = xy[i1]
            x2, y2 = xy[i2]
            lo_x = max(int(math.floor(min(x0, x1, x2))), 0)
            hi_x = min(int(math.ceil(max(x0, x1, x2))) + 1, S)
            lo_y = max(int(math.floor(min(y0, y1, y2))), 0)
            hi_y = min(int(math.ceil(max(y0, y1, y2))) + 1, S)
            if lo_x >= hi_x or lo_y >= hi_y:
                continue
            den = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2)
            if abs(den) < 1e-12:
                continue
            gx, gy = np.meshgrid(np.arange(lo_x, hi_x) + 0.5,
                                 np.arange(lo_y, hi_y) + 0.5)
            w0 = ((y1 - y2) * (gx - x2) + (x2 - x1) * (gy - y2)) / den
            w1 = ((y2 - y0) * (gx - x2) + (x0 - x2) * (gy - y2)) / den
            w2 = 1.0 - w0 - w1
            inside = (w0 >= 0) & (w1 >= 0) & (w2 >= 0)
            if not inside.any():
                continue
            depth = w0 * z[i0] + w1 * z[i1] + w2 * z[i2]
            sub = zbuf[lo_y:hi_y, lo_x:hi_x]
            win = inside & (depth > sub)
            if not win.any():
                continue
            sub[win] = depth[win]
            rgb[lo_y:hi_y, lo_x:hi_x][win] = tone(n[i], tones)

        owner = {}
        for i in range(len(f)):
            for k in range(3):
                u, w = int(f[i, k]), int(f[i, (k + 1) % 3])
                owner.setdefault((u, w) if u < w else (w, u), []).append(i)
        cos_lim = math.cos(math.radians(feature_deg))
        for (u, w), tris in owner.items():
            if len(tris) != 2:
                continue
            i, j = tris
            if n[i, 2] <= 0.0 and n[j, 2] <= 0.0:
                continue
            if float(np.dot(n[i], n[j])) < cos_lim:
                creases.append((tuple(xy[u]), tuple(xy[w])))

    shaded = Image.fromarray(rgb, 'RGBA')
    d = ImageDraw.Draw(shaded)
    for seg in creases:
        d.line(list(seg), fill=INK, width=max(1, int(round(ew * 0.55))))

    # The outline is the shape's own silhouette DILATED, not a stroke along the
    # polygon boundary: a stroke is centred on the edge, so half of it falls
    # inside and the faces painted next cover it. Growing the alpha mask by the
    # pen width and filling that with ink puts the whole outline outside, and it
    # follows the union of both solids with no gaps at the joints.
    mask = shaded.split()[3].point(lambda a: 255 if a > 96 else 0)
    grown = mask
    k = 2 * ew + 1
    while k > 1:
        step = min(k, 9)                 # MaxFilter is slow with a wide kernel
        grown = grown.filter(ImageFilter.MaxFilter(step if step % 2 else step - 1))
        k -= (step if step % 2 else step - 1) - 1
    out = Image.new('RGBA', (S, S), (0, 0, 0, 0))
    out.paste(Image.new('RGBA', (S, S), INK), (0, 0), grown)
    out.paste(shaded, (0, 0), shaded)
    return out.resize((size, size), Image.LANCZOS)


def main():
    os.makedirs(OUT, exist_ok=True)
    sizes = ((16, 0.075), (32, 0.060), (64, 0.050))
    for s, ew in sizes:
        render(s, edge_frac=ew).save(os.path.join(OUT, '%dx%d.png' % (s, s)))
        print('wrote %dx%d.png' % (s, s))
    prev = os.path.join(HERE, '_preview')
    os.makedirs(prev, exist_ok=True)
    for s, ew in sizes:
        render(s, edge_frac=ew).resize((256, 256), Image.NEAREST).save(
            os.path.join(prev, 'wrap_%d_zoom.png' % s))
    print('previews in %s' % prev)


if __name__ == '__main__':
    main()
