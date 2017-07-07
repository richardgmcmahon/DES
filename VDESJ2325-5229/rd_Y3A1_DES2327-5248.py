from __future__ import (print_function, division)

# standard library
import os
import sys
import time


# 3rd party
import numpy as np
from matplotlib import pyplot as plt

from astropy.table import Table

# private



if __name__ == "__main__":

    import argparse
    import configparser

    # inpath = './'
    inpath = '/data/desardata3/Y3A1/r2587/DES2327-5248/p01/cat/'

    filename_des = 'DES2327-5248_r2587p01.fits'

    filename_desxmatch = 'DES2327-5248_r2587p01_matched.fits'

    infile_des = inpath + filename_des
    des = Table.read(infile_des)
    des.info()

    des.info('stats')

    infile_desxmatch = inpath + filename_desxmatch
    desxmatch = Table.read(infile_desxmatch)
    desxmatch.info()
    desxmatch.info('stats')
