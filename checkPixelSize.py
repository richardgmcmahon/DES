from __future__ import (print_function, division)
"""
todo check x, y orientation of images versus numpy orientation
using an image that has unequal axes e.g. VHS tile

History:

2016-04-01: rgm: Original


"""


import os
import sys
import time

import numpy as np

from astropy.io import fits
from astropy import wcs


def wcs_PixelScale(AstWCS, x, y, debug=False):
    """
    simple determination of the pixel scale by computing the change in
    RA, Dec by one pixel centered on a specific pixel

         4
      1  0  2
         3

    """

    ra1, dec1 = AstWCS.wcs_pix2world(x - 0.5, y, 1)
    ra2, dec2 = AstWCS.wcs_pix2world(x + 0.5, y, 1)

    debug = False
    if debug:
        print('dec1, dec2:', dec1, dec2)
    dec = (dec1 + dec2) / 2.0
    if debug:
        print('Dec, cos(Dec):', dec, np.cos(np.deg2rad(dec)))
    RAScale = (ra1 - ra2) * 3600.0 * np.cos(np.deg2rad(dec))

    if debug:
        print('RAScale:', RAScale)

    ra3, dec3 = AstWCS.wcs_pix2world(x, y - 0.5, 1)
    ra4, dec4 = AstWCS.wcs_pix2world(x, y + 0.5, 1)

    DecScale = (dec4 - dec3) * 3600.0

    print('x, y, RAScale, DecScale, AspectRatio:',
          x, y, RAScale, DecScale, RAScale / DecScale)

    return RAScale, DecScale


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument("--file",
                        help="input file")

    parser.add_argument("--ext", type=int, default=1,
                        help="image extension 0 to n; [Default: 1]")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    infile = args.file
    ext = args.ext

    hdulist = fits.open(infile)
    hdulist.info()

    AstWCS = wcs.WCS(hdulist[ext].header)

    NAXIS1 = hdulist[ext].header['NAXIS1']
    NAXIS2 = hdulist[ext].header['NAXIS2']
    print('NAXIS1:', NAXIS1)
    print('NAXIS2:', NAXIS2)

    print('Processing:', infile)
    x = (NAXIS1 + 1) / 2.0
    y = (NAXIS2 + 1) / 2.0
    wcs_PixelScale(AstWCS, x, y)
    wcs_PixelScale(AstWCS, 1, 1)
    wcs_PixelScale(AstWCS, 1, NAXIS2)
    wcs_PixelScale(AstWCS, NAXIS1, 1)
    wcs_PixelScale(AstWCS, NAXIS1, NAXIS2)
