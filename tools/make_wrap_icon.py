"""Build the Wrap button's icon set from the artwork in tools/artwork.

    py -3.10 tools/make_wrap_icon.py

THE TWO SOURCE IMAGES
    Wrap_Icon_Fusion_Render.png is the icon: a blue offset shell with the body it
    came from in front of it, and a dark strip between the two.

    Icono_Wrap_Verde.png is the same render with that strip painted GREEN. It is
    a colour key rather than a description, so which pixels the strip occupies
    needs no guessing here. Those pixels are removed.

WHAT FUSION'S OWN ICONS ACTUALLY DO
    Measured off the installed set (Fusion/UI/FusionUI/Resources/solid/shell,
    solid/thicken, Modeling/OffsetFaces), not guessed:

    - Each command ships SIX pngs plus HiDPI twins: 16x16 / 32x32 / 64x64 for the
      light theme and the same names with "-dark" for the dark one, plus "@2x"
      copies at double resolution. A folder holding only the light three leaves
      Fusion nothing to load on the dark theme.
    - The two themes are not one image reused. The LIGHT variant is drawn in
      #919191..#F3F3F3 greys and #3F85B7..#89CBFA blues and carries a 2 px outline
      (at 64) in #666666. The DARK variant is lighter throughout -
      #BBBBBB..#F3F3F3 and #67B1E6..#B2DFFF - and carries NO outline at all: on a
      dark ribbon the silhouette does that work by itself.
    - They fill the canvas. shell and OffsetFaces span all 64 px of theirs; the
      margin an icon can afford is a pixel or two, not a tenth of the width.

    The version this replaces used one image for both themes, a 9 % margin and an
    outline six times too thick in the ribbon's own colour. That is what looked
    like low resolution.
"""
import os

import numpy as np
from PIL import Image
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
ART = os.path.join(HERE, 'artwork')
OUT = os.path.join(os.path.dirname(HERE), 'resources', 'Wrap')

RENDER = os.path.join(ART, 'Wrap_Icon_Fusion_Render.png')
GREEN = os.path.join(ART, 'Icono_Wrap_Verde.png')

# darkest and lightest of each material, sampled from the icons named above
THEMES = {
    'light': {'grey': ((0x91, 0x91, 0x91), (0xF3, 0xF3, 0xF3)),
              'blue': ((0x3F, 0x85, 0xB7), (0x89, 0xCB, 0xFA)),
              'ink': (0x66, 0x66, 0x66),
              'outline_px_at_64': 2.0},
    'dark': {'grey': ((0xBB, 0xBB, 0xBB), (0xF3, 0xF3, 0xF3)),
             'blue': ((0x67, 0xB1, 0xE6), (0xB2, 0xDF, 0xFF)),
             'ink': None,
             'outline_px_at_64': 0.0},
}

SUPERSAMPLE = 8
SPECK = 40                 # px at source scale: antialiasing crumbs, not shape
MARGIN = 0.02              # of the canvas, per side


def artwork():
    """The render with the green-keyed strip removed, cleaned and cropped."""
    a = np.array(Image.open(RENDER).convert('RGBA'))
    g = np.array(Image.open(GREEN).convert('RGBA')).astype(int)
    if a.shape != g.shape:
        raise RuntimeError('the two source images are different sizes')
    green = ((g[:, :, 1] > g[:, :, 0] + 40) & (g[:, :, 1] > g[:, :, 2] + 40)
             & (g[:, :, 3] > 128))
    art = a.copy()
    art[green] = (0, 0, 0, 0)

    # removing the strip leaves a couple of dozen four-pixel crumbs along its
    # antialiased border; each would survive the downsample as a speck
    mask = art[:, :, 3] > 128
    lbl, n = ndimage.label(mask)
    if n > 1:
        areas = ndimage.sum(mask, lbl, range(1, n + 1))
        keep = np.zeros(n + 1, bool)
        for i, area in enumerate(areas, start=1):
            keep[i] = area >= SPECK
        art[~keep[lbl]] = (0, 0, 0, 0)

    op = art[:, :, 3] > 40
    ys, xs = np.nonzero(op)
    return art[ys.min():ys.max() + 1, xs.min():xs.max() + 1]


def recolour(rgba, theme):
    """Move the render onto Fusion's palette while keeping its own shading.

    Each material is remapped by LUMINANCE rather than replaced with a flat
    colour: the darkest pixel of the blue becomes the palette's darkest blue, the
    lightest its lightest, and everything between rides along. The render's own
    modelling survives; only the range it occupies moves. Iago's blue sits at
    #618AAA and his grey at #767778, both far darker than Fusion's, so without
    this the icon reads as a different family however right its shape.
    """
    pal = THEMES[theme]
    out = rgba.copy()
    op = rgba[:, :, 3] > 0
    px = rgba[:, :, :3].astype(float)
    blue = op & (px[:, :, 2] > px[:, :, 0] + 20)
    grey = op & ~blue

    for sel, key in ((grey, 'grey'), (blue, 'blue')):
        if not sel.any():
            continue
        lum = px[sel].mean(axis=1)
        lo, hi = float(lum.min()), float(lum.max())
        t = (lum - lo) / max(hi - lo, 1e-6)
        dark, light = (np.array(c, float) for c in pal[key])
        out[sel, :3] = np.clip(dark + (light - dark) * t[:, None],
                               0, 255).astype(np.uint8)
    return out


def build(size, theme, scale=1):
    """One png. `scale` 2 makes the @2x twin, at double resolution."""
    px = size * scale
    S = px * SUPERSAMPLE
    pal = THEMES[theme]

    art = recolour(artwork(), theme)
    inner = int(round(S * (1.0 - 2.0 * MARGIN)))
    small = Image.fromarray(art).resize((inner, inner), Image.LANCZOS)
    canvas = Image.new('RGBA', (S, S), (0, 0, 0, 0))
    off = (S - inner) // 2
    canvas.paste(small, (off, off), small)
    rgba = np.array(canvas)
    mask = rgba[:, :, 3] > 128

    # the outline is a DILATION, not a stroke: a stroke is centred on the edge,
    # so half of it lands inside the artwork and gets covered
    pen = pal['outline_px_at_64'] * px / 64.0 * SUPERSAMPLE
    if pal['ink'] is not None and pen >= 1.0:
        grown = ndimage.distance_transform_edt(~mask) <= pen
        out = np.zeros_like(rgba)
        out[grown] = tuple(pal['ink']) + (255,)
        out[mask] = rgba[mask]
    else:
        out = rgba

    return Image.fromarray(out).resize((px, px), Image.LANCZOS)


def main():
    os.makedirs(OUT, exist_ok=True)
    written = []
    for size in (16, 32, 64):
        for theme in ('light', 'dark'):
            suffix = '' if theme == 'light' else '-dark'
            for scale in (1, 2):
                if scale == 2 and size == 64:
                    continue           # Fusion ships no 64x64@2x either
                at = '@2x' if scale == 2 else ''
                name = '%dx%d%s%s.png' % (size, size, suffix, at)
                build(size, theme, scale).save(os.path.join(OUT, name))
                written.append(name)
    print('wrote %d files to %s' % (len(written), OUT))
    for n in sorted(written):
        print('   %s' % n)

    prev = os.path.join(HERE, '_preview')
    os.makedirs(prev, exist_ok=True)
    for theme in ('light', 'dark'):
        for size in (16, 32, 64):
            build(size, theme).resize((256, 256), Image.NEAREST).save(
                os.path.join(prev, 'wrap_%s_%d.png' % (theme, size)))


if __name__ == '__main__':
    main()
