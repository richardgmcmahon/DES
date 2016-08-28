from __future__ import print_function, division

import sys

import numpy as np

import matplotlib.pyplot as plt

from astropy.table import Table

sys.path.append('/home/rgm/soft/python/lib/')
from librgm import plotid


def flags_stats(table):
    """


    """

    print()
    print('Number of data values:', len(table))

    for band in ['g', 'r', 'i', 'z', 'Y']:

        FLAG = 'DES ' + band + 'Flag'

        itest = (table['DES ' + band + 'Flag'] == 0)
        print('Number zero = 0 FLAGS_' + band + ': ',
              len(table[itest]))

        itest = (table['DES ' + band + 'Flag'] > 0)
        print('Number zero > 0 FLAGS_' + band + ': ',
              len(table[itest]))

        itest = (table['DES ' + band + 'Flag'] < 0)
        print('Number zero < 0 FLAGS_' + band + ': ',
              len(table[itest]))

        print('Minimum and Maximum FLAGS' + band + ': ',
              np.min(table[FLAG]), np.max(table[FLAG]))

        print()

    itest = ((table['DES gFlag'] == 0) & (table['DES rFlag'] == 0) &
             (table['DES iFlag'] == 0) & (table['DES zFlag'] == 0) &
             (table['DES YFlag'] == 0))

    print('Number FLAG = 0  bands GRIZY: ',
          len(table[itest]))

    itest = ((table['DES gFlag'] <= 3) & (table['DES rFlag'] <= 3) &
             (table['DES iFlag'] <= 3) & (table['DES zFlag'] <= 3) &
             (table['DES YFlag'] <= 3))

    print('Number FLAG <= 3  ALL bands GRIZY: ',
          len(table[itest]))

    itest = ((table['DES gFlag'] > 0) & (table['DES rFlag'] > 0) &
             (table['DES iFlag'] > 0) & (table['DES zFlag'] > 0) &
             (table['DES YFlag'] > 0))

    print('Number zero > 0  ALL bands GRIZY: ',
          len(table[itest]))

    itest = ((table['DES gFlag'] < 0) & (table['DES rFlag'] < 0) &
             (table['DES iFlag'] < 0) & (table['DES zFlag'] < 0) &
             (table['DES YFlag'] < 0))

    print('Number zero < 0  ALL bands GRIZY: ',
          len(table[itest]))

    itest = ((table['DES gFlag'] < 0) & (table['DES rFlag'] < 0) &
             (table['DES iFlag'] < 0) & (table['DES zFlag'] < 0) &
             (table['DES YFlag'] < 0))

    print('Number zero < 0  ALL bands GRIZY: ',
          len(table[itest]))


    itest = ((table['DES gFlag'] = 0) | (table['DES rFlag'] = 0) |
             (table['DES iFlag'] = 0) | (table['DES zFlag'] = 0) |
             (table['DES YFlag'] = 0))

    print('Number FLAGS_GRIZY = 0 in at least on band: ',
          len(table[itest]))

    FLAG = 'DES gFlag'
    xdata = table[FLAG]
    print(np.min(xdata), np.max(xdata))
    xrange = np.max(xdata) - np.min(xdata)
    print(xrange)
    bins = 7
    range = [-0.5, bins - 0.5]
    ndata = len(xdata)
    plt.hist(xdata-0.5, bins=bins, range=range,
             align='mid', histtype='step',
             label=FLAG + ': ' + str(ndata))

    FLAG = 'DES rFlag'
    xdata = table[FLAG]
    print(np.min(xdata), np.max(xdata))
    xrange = np.max(xdata) - np.min(xdata)
    print(xrange)
    bins = 7
    ndata = len(xdata)
    plt.hist(xdata-0.5, bins=bins, range=range,
             align='mid', histtype='step',
             label=FLAG + ': ' + str(ndata))


    FLAG = 'DES iFlag'
    xdata = table[FLAG]
    print(np.min(xdata), np.max(xdata))
    xrange = np.max(xdata) - np.min(xdata)
    print(xrange)
    bins = 7
    ndata = len(xdata)
    plt.hist(xdata-0.5, bins=bins, range=range,
             align='mid', histtype='step',
             label=FLAG + ': ' + str(ndata))


    FLAG = 'DES zFlag'
    xdata = table[FLAG]
    print(np.min(xdata), np.max(xdata))
    xrange = np.max(xdata) - np.min(xdata)
    print(xrange)
    bins = 7
    ndata = len(xdata)
    plt.hist(xdata-0.5, bins=bins, range=range,
             align='mid', histtype='step',
             label=FLAG + ': ' + str(ndata))


    FLAG = 'DES YFlag'
    xdata = table[FLAG]
    print(np.min(xdata), np.max(xdata))
    xrange = np.max(xdata) - np.min(xdata)
    print(xrange)
    bins = 7
    ndata = len(xdata)
    plt.hist(xdata-0.5, bins=bins, range=range,
             align='mid', histtype='step',
             label=FLAG + ': ' + str(ndata))


    # plt.hist(xdata, bins, histtype='step')
    # plt.hist(xdata, histtype='step')
    # plt.bar(xdata-0.5, 4, width=0.8)
    # plt.xticks(range(4))

    print(range)
    plt.xlim(range)
    plt.xlabel('FLAGS_GRIZY')
    plt.ylabel('Number')
    plt.title(filename)
    plt.legend(fontsize="medium")
    plotid.plotid()
    plt.show()

    return


if __name__ == "__main__":


    inpath = "/home/rgm/Projects/DES/OzDES/"
    filename = "OZDES_QSO_20160627.fits"
    infile = inpath + filename

    table = Table.read(infile)
    table.info()
    table.info('stats')

    flags_stats(table)

    xsize = 8
    ysize = 8
    plt.figure(figsize=(xsize, ysize))

    xdata = table['REDSHIFT']
    ydata = table['DES coadd i']

    plt.scatter(xdata, ydata)

    plt.show()

    plt.figure(figsize=(xsize, ysize))

    itest = (table['DES coadd i'] < 90) & (table['DES coadd z'] < 90)
    xdata = table['REDSHIFT'][itest]
    ydata = table['DES coadd i'][itest] - table['DES coadd z'][itest]

    ndata = len(xdata)
    plt.scatter(xdata, ydata, color='blue', edgecolor='none',
                label=str(ndata))
    plt.title(filename)
    plt.xlabel('Redshift')
    plt.ylabel('i - z [DES coadd]')
    plt.xlim(0.0, 5.0)
    plt.legend()

    print(np.median(xdata), np.median(ydata), len(xdata), len(ydata))

    plt.show()


    itest = (table['DES coadd g'] < 90) & (table['DES coadd i'] < 90)
    xdata = table['REDSHIFT'][itest]
    ydata = table['DES coadd g'][itest] - table['DES coadd i'][itest]

    plt.scatter(xdata, ydata)

    plt.show()

    itest = (table['W1'] < 90) & (table['W2'] < 90)
    xdata = table['REDSHIFT'][itest]
    ydata = table['W1'][itest] - table['W2'][itest]

    plt.scatter(xdata, ydata)

    plt.show()


    itest = (table['W1'] < 90) & (table['W2'] < 90)
    xdata = table['REDSHIFT'][itest]
    ydata = table['W1'][itest] - table['W2'][itest]

    plt.scatter(xdata, ydata)

    plt.show()
