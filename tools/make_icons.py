"""Generate the toolbar icon set for the Mesh to Part add-in.

Left half of the silhouette is drawn as faceted mesh, right half as a smooth
solid. Rendered at 8x and downsampled so it stays clean at 16 px.
"""
import math, os
from PIL import Image, ImageDraw

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "MeshToPart")
os.makedirs(OUT, exist_ok=True)

MESH_FILL = (208, 214, 221, 255)
MESH_EDGE = (94, 106, 120, 255)
SOLID_FILL = (44, 122, 168, 255)
SOLID_EDGE = (24, 78, 112, 255)


def silhouette(cx, cy, r, n=96):
    pts = []
    for i in range(n):
        t = 2.0 * math.pi * i / n
        rad = r * (1.0 + 0.13 * math.sin(t * 2.0 + 0.6) - 0.07 * math.cos(t * 3.0))
        pts.append((cx + rad * math.cos(t), cy + rad * math.sin(t) * 0.92))
    return pts


def render(size, facets, lwf=0.018, rf=0.40):
    S = size * 8
    img = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)
    cx = cy = S / 2.0
    r = S * rf
    poly = silhouette(cx, cy, r)

    left = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    dl = ImageDraw.Draw(left)
    dl.polygon(poly, fill=MESH_FILL, outline=MESH_EDGE, width=max(1, int(S * lwf)))
    lw = max(1, int(S * lwf))
    # facet chords across the whole blob, later clipped to the left half
    ring = silhouette(cx, cy, r * 0.52, facets)
    outer = silhouette(cx, cy, r, facets)
    for i in range(facets):
        dl.line([outer[i], ring[i]], fill=MESH_EDGE, width=lw)
        dl.line([ring[i], ring[(i + 1) % facets]], fill=MESH_EDGE, width=lw)
        dl.line([outer[i], outer[(i + 1) % facets]], fill=MESH_EDGE, width=lw)
    mask_l = Image.new("L", (S, S), 0)
    ImageDraw.Draw(mask_l).rectangle([0, 0, int(cx) - 1, S], fill=255)
    img.paste(left, (0, 0), Image.composite(left.split()[3], Image.new("L", (S, S), 0), mask_l))

    right = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    dr = ImageDraw.Draw(right)
    dr.polygon(poly, fill=SOLID_FILL, outline=SOLID_EDGE, width=lw)
    mask_r = Image.new("L", (S, S), 0)
    ImageDraw.Draw(mask_r).rectangle([int(cx), 0, S, S], fill=255)
    img.paste(right, (0, 0), Image.composite(right.split()[3], Image.new("L", (S, S), 0), mask_r))

    # the seam
    ys = [p[1] for p in poly if abs(p[0] - cx) < r * 0.35]
    if ys:
        d.line([(cx, min(ys)), (cx, max(ys))], fill=SOLID_EDGE, width=lw)

    return img.resize((size, size), Image.LANCZOS)


for s, f, w, rf in ((16, 5, 0.048, 0.45), (32, 7, 0.026, 0.43), (64, 9, 0.018, 0.41)):
    render(s, f, w, rf).save(os.path.join(OUT, "%dx%d.png" % (s, s)))
    print("wrote %dx%d.png" % (s, s))
