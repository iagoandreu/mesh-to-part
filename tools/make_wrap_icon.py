"""Generate the toolbar icon for the Wrap (offset) command.

    py -3.10 tools/make_wrap_icon.py

The idea is 3-matic's Wrap icon - a skin pulled over a block, its hem hanging
loose and rolled at one corner - drawn in the style the rest of this add-in
already uses: light grey body, grey-blue edges, the house blue as the accent,
thin strokes rather than 3-matic's heavy black outline.

The two buttons shared one icon before this, which is the worst case for a
toolbar: the only way to tell Convert from Wrap was to read the label.

Rendered at 8x and downsampled, same as make_icons.py. What survives 16 px is the
silhouette, so the hem roll deliberately projects PAST the box on the left and the
hem flips up past it on the right - that is what makes this button a different
shape from the round Convert one at a glance, rather than a different colour.
"""
import math
import os

from PIL import Image, ImageDraw

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   'resources', 'Wrap')

# The body and skin colours are the ones make_icons.py already uses. The EDGE is
# not: at 1.8 % of the canvas in a low-contrast blue, every band merged into one
# mass and the icon read as a hat. Profiling the 3-matic original shows it
# separates each band with a stroke about 4.7 % of the drawing height, and that
# separation is what makes the shape readable at all. This keeps the separation
# and softens the colour from pure black to a dark navy, which is the Fusion half
# of the brief.
BODY_TOP = (226, 231, 236, 255)
BODY_LEFT = (208, 214, 221, 255)
BODY_RIGHT = (185, 193, 203, 255)

SKIN_TOP = (96, 166, 212, 255)
SKIN_LEFT = (52, 130, 178, 255)
SKIN_RIGHT = (38, 104, 146, 255)

INK = (27, 42, 56, 255)


def lerp(p, q, t):
    return (p[0] + (q[0] - p[0]) * t, p[1] + (q[1] - p[1]) * t)


def bezier(pts, n=64):
    """Sample a de Casteljau curve through the control points."""
    out = []
    for i in range(n + 1):
        t = i / float(n)
        cur = list(pts)
        while len(cur) > 1:
            cur = [lerp(cur[k], cur[k + 1], t) for k in range(len(cur) - 1)]
        out.append(cur[0])
    return out


def tube(d, path, radii, _unused, fill, edge, w):
    """A rolled hem: a stack of circles, one radius per path point.

    A constant-width stroke ends in a round cap, and a round cap in the middle of
    the box reads as a bar laid across it. Per-point radii let the hem thin out
    towards the far corner the way the original's does.

    Each radius is floored at 1.2 w. Below that the outline meets itself and the
    thin end goes solid dark - which is how an early curl came out as a 27 px
    black bar in the profile.

    Two passes, edge then fill: stroking the whole path twice cannot leave a seam
    where segments meet, which outlining each segment separately does.
    """
    for pas in (0, 1):
        col = edge if pas == 0 else fill
        for c, rr in zip(path, radii):
            rr = max(rr, 1.2 * w)
            if pas == 0:
                rr += w
            d.ellipse([c[0] - rr, c[1] - rr, c[0] + rr, c[1] + rr], fill=col)


def icon(size, lwf=0.041, simple=False):
    """Proportions and stroke weight measured off the 3-matic icon, not guessed.

    Profiling a column down the middle of it: the top face is 27 % of the drawing
    height and the front face 67 %, and that front face divides 42 % skin / 28 %
    rolled hem / 30 % bare block. Every band is separated from the next by a
    stroke about 4.7 % of the height - which is thick, and is the whole reason
    the original stays legible at 16 px.
    """
    S = size * 8
    img = Image.new('RGBA', (S, S), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)
    w = max(1, int(round(S * lwf)))

    cx = S * 0.475
    a = S * 0.310                       # half width of the top rhombus
    b = S * 0.100                       # half height of it
    hh = S * 0.660                      # front face height
    ty = S * 0.048 + b

    L = (cx - a, ty)
    T = (cx, ty - b)
    R = (cx + a, ty)
    B = (cx, ty + b)
    Lb = (L[0], L[1] + hh)
    Bb = (B[0], B[1] + hh)
    Rb = (R[0], R[1] + hh)

    r = hh * 0.092                      # roll radius: 2r is 18 % of the face
    # 47 % skin / 24 % roll / 29 % bare block is what a column down the middle of
    # the original gives. Earlier drafts put the roll at 42 % with 40 % of block
    # below it, which sat the roll at mid-height and made a tray for the box to
    # stand on instead of a hem near its foot.
    band = 0.50                         # of hh, skin above the roll

    # ---- the block ------------------------------------------------------- #
    d.polygon([L, B, Bb, Lb], fill=BODY_LEFT)
    d.polygon([B, R, Rb, Bb], fill=BODY_RIGHT)

    # ---- the skin, down to where the roll starts ------------------------- #
    hlL = (L[0], L[1] + hh * band)
    hlB = (B[0], B[1] + hh * (band + 0.05))
    hrR = (R[0], R[1] + hh * (band - 0.12))
    d.polygon([L, B, hlB, hlL], fill=SKIN_LEFT)
    d.polygon([B, R, hrR, hlB], fill=SKIN_RIGHT)
    d.polygon([L, T, R, B], fill=SKIN_TOP)

    # ---- every edge of the block, BEFORE the roll ------------------------ #
    # The roll passes in FRONT of the box where it lifts, so the box's own edges
    # have to go down first and let the roll cover them. Drawing them afterwards
    # cut the roll in two at the right corner and it read as a belt passing
    # behind the box.
    # At 16 px every stroke is under a pixel and they bleed into each other, so
    # the internal vertical edge goes: the three faces still read from their
    # fills, and losing one line is cheaper than losing all of them to mud.
    edges = [[L, Lb], [Lb, Bb], [Bb, Rb], [R, Rb], [L, T], [T, R], [L, B], [B, R]]
    if not simple:
        edges.append([B, Bb])
    for seg in edges:
        d.line(seg, fill=INK, width=w)

    # ---- the hem: ONE roll, an offset of the hem across both faces -------- #
    # Profiling the original column by column settles what this is. At 25 % and
    # 45 % across it reads skin / line / roll / line / block; at 60 % the line
    # between skin and roll has vanished and they are one 29 px band; at 85 % the
    # line is back with a THINNER roll below it. So it is a single hem running
    # the whole way round, tangent to the skin, tapering to the right - not a bar
    # laid across the box, and not two separate blobs.
    #
    # Building it as an offset of the hem line rather than from hand-placed
    # bezier controls is what stopped it drifting into a ribbon or a hat brim.
    # How far the roll pokes past the box. Measured on the original it is about
    # 9 % of the box width per side; at 17 % - where an earlier draft had it - the
    # roll reads as a collar and the box looks like it is standing on a tray.
    ov = S * 0.004
    sl_l = (hlB[1] - hlL[1]) / (hlB[0] - hlL[0])
    sl_r = (hrR[1] - hlB[1]) / (hrR[0] - hlB[0])
    p_l = (hlL[0] - ov, hlL[1] - sl_l * ov + r)
    p_m = (hlB[0], hlB[1] + r)
    p_r = (hrR[0] + ov, hrR[1] + sl_r * ov + r * 0.62)

    roll = [lerp(p_l, p_m, i / 44.0) for i in range(45)]
    roll += [lerp(p_m, p_r, i / 44.0) for i in range(1, 45)]
    # Fat on the left face, thinning away to almost nothing on the right. Holding
    # the width right round made a ledge the box appeared to stand on; the
    # original rolls its hem at the near-left corner and lets the far side lie
    # flat, which is also what a sheet actually does when you pull it over a block.
    radii = [r] * 30
    radii += [r + (r * 0.18 - r) * (i / 58.0) for i in range(59)]
    tube(d, roll, radii, None, SKIN_LEFT, INK, w)

    # the far corner lifted off the block: short, low, just past the silhouette
    if not simple:
        flap = bezier([(R[0] - a * 0.38, hrR[1] + r * 0.30),
                       (R[0] - a * 0.06, hrR[1] + r * 0.05),
                       (R[0] + S * 0.018, hrR[1] - hh * 0.10)], 40)
        tube(d, flap, [r * 0.62] * 41, None, SKIN_RIGHT, INK, w)

    d.line([hlL, hlB, hrR], fill=INK, width=w, joint='curve')

    return img.resize((size, size), Image.LANCZOS)


def main():
    os.makedirs(OUT, exist_ok=True)
    for s, lw, simple in ((16, 0.040, True), (32, 0.042, False),
                          (64, 0.041, False)):
        icon(s, lw, simple).save(os.path.join(OUT, '%dx%d.png' % (s, s)))
        print('wrote %dx%d.png' % (s, s))
    # magnified copies to judge by while tuning; kept OUT of resources/ so they
    # are not shipped into Fusion's AddIns folder by install.ps1
    prev = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_preview')
    os.makedirs(prev, exist_ok=True)
    icon(16, 0.040, True).resize((256, 256), Image.NEAREST).save(
        os.path.join(prev, 'wrap_16_zoom.png'))
    icon(64, 0.041).resize((256, 256), Image.NEAREST).save(
        os.path.join(prev, 'wrap_64_zoom.png'))
    print('previews written to %s' % prev)


if __name__ == '__main__':
    main()
