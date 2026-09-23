"""Build the Wrap button's icon set from the artwork in tools/artwork.

    py -3.10 tools/make_icons.py

THE TWO SOURCE IMAGES
    Wrap_Icon_Fusion_Render.png is the icon: a blue offset shell with the body it
    came from in front of it, and a dark strip between the two.

    Icono_Wrap_Verde.png is the same render with that strip painted GREEN. It is
    there to say, without any guessing on this side, exactly which pixels the
    strip occupies - a colour key rather than a description. Those pixels are
    removed from the render, which leaves a gap and, with it, two separate pieces.

THE OUTLINE IS A DILATION, NOT A STROKE
    The shape is then expanded pixel by pixel - a euclidean distance transform, so
    the growth is a disc in every direction rather than a square - and everything
    the expansion adds is filled with the ribbon's own #3B4453.

    That does three things at once: it closes the gap the green left, so the icon
    is a single connected piece again rather than two; it gives every contour a
    thick dark edge in exactly the shape of the icon; and because the colour is
    the ribbon's, the outline disappears on Fusion's dark theme and only shows
    where it is needed, on a light one.

    A stroke along the silhouette would not do this. A stroke is centred on the
    edge, so half of it lands inside and gets covered by the artwork, and it has
    no idea the two pieces should be joined.
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

INK = (59, 68, 83, 255)            # #3B4453, the Fusion ribbon's own dark grey

# outline thickness as a fraction of the icon, per size. Smaller icons need
# proportionally more or the outline lands under a pixel and vanishes.
THICKNESS = {16: 0.085, 32: 0.070, 64: 0.060}

SUPERSAMPLE = 8
SPECK = 40                          # px at working scale; antialiasing crumbs


def artwork():
    """The render with the green-keyed strip removed, and the crumbs cleaned."""
    a = np.array(Image.open(RENDER).convert('RGBA'))
    g = np.array(Image.open(GREEN).convert('RGBA')).astype(int)
    green = ((g[:, :, 1] > g[:, :, 0] + 40) & (g[:, :, 1] > g[:, :, 2] + 40)
             & (g[:, :, 3] > 128))
    if a.shape != g.shape:
        raise RuntimeError('the two source images are different sizes')
    art = a.copy()
    art[green] = (0, 0, 0, 0)
    return art, int(green.sum())


def despeckle(alpha, min_area):
    """Drop islands smaller than min_area.

    Removing the strip leaves a couple of dozen four-pixel crumbs along its
    antialiased border. Each one would grow into a dark dot of its own under the
    dilation, so they go first.
    """
    mask = alpha > 128
    lbl, n = ndimage.label(mask)
    if n <= 1:
        return mask
    keep = np.zeros(n + 1, bool)
    for i, area in enumerate(ndimage.sum(mask, lbl, range(1, n + 1)), start=1):
        keep[i] = area >= min_area
    return keep[lbl]


def render(size):
    art, green_px = artwork()
    S = size * SUPERSAMPLE
    small = Image.fromarray(art).resize((S, S), Image.LANCZOS)
    rgba = np.array(small)

    mask = despeckle(rgba[:, :, 3], SPECK)
    rgba[~mask] = (0, 0, 0, 0)

    # the pen in working pixels: a fraction of the FINAL icon, scaled up
    pen = max(1, int(round(size * THICKNESS[size] * SUPERSAMPLE)))

    # the expansion itself: every pixel within `pen` of the shape, measured as a
    # true distance so the growth is round rather than square
    dist = ndimage.distance_transform_edt(~mask)
    grown = dist <= pen

    out = np.zeros((S, S, 4), np.uint8)
    out[grown] = INK
    out[mask] = rgba[mask]

    lbl, n = ndimage.label(grown)
    img = Image.fromarray(out).resize((size, size), Image.LANCZOS)
    return img, {'green_px': green_px, 'pen_px': pen, 'pieces': n}


def main():
    os.makedirs(OUT, exist_ok=True)
    prev = os.path.join(HERE, '_preview')
    os.makedirs(prev, exist_ok=True)
    for size in (16, 32, 64):
        img, info = render(size)
        img.save(os.path.join(OUT, '%dx%d.png' % (size, size)))
        img.resize((256, 256), Image.NEAREST).save(
            os.path.join(prev, 'wrap_%d_zoom.png' % size))
        print('%2d px: pen %3d px at %dx supersample -> %d connected piece%s'
              % (size, info['pen_px'], SUPERSAMPLE, info['pieces'],
                 '' if info['pieces'] == 1 else 's  *** should be 1'))


if __name__ == '__main__':
    main()
