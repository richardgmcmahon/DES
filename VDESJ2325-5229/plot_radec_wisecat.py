import sys

sys.path.append('/home/rgm/soft/python/lib/')
from librgm.plotid import plotid



def plot_radec_wisecat(data=None,
                       source=None,
                       radec_centre=None,
                       xrange = [-2.5, 2.5],
                       yrange = [-2.5, 2.5],
                       overplot=True,
                       plt=None,
                       showplot=False,
                       debug=False):
    """


    """

    import numpy as np
    from matplotlib.patches import Circle

    print('filename:', data.meta['filename'])

    ra0 = radec_centre[0]
    dec0 = radec_centre[1]

    print('Number of rows:', len(data))
    if debug:
        data.info()
        data.info('stats')
        print('zoom:', zoom)
        print('filename:', data.meta['filename'])
        print('ra0', ra0)
        print('dec0', dec0)
        print('xrange:', xrange)
        print('yrange:', yrange)

    ra = data['RAJ2000']
    dec = data['DEJ2000']

    ra_min = np.min(ra)
    ra_max = np.max(ra)
    delta_ra = (ra - ra0) * 3600.0 * np.cos(np.deg2rad(dec0))

    if debug:
        print('ndata', len(ra))
        print('ra_min:', ra_min)
        print('ra_max:', ra_max)
        print('Delta RA range:',
              np.min(delta_ra), np.max(delta_ra))

    dec_min = np.min(dec)
    dec_max = np.max(dec)
    delta_dec = (dec - dec0) * 3600.0

    if debug:
        print('dec_min:', dec_min)
        print('dec_max:', dec_max)
        print('Delta Dec range:',
              np.min(delta_dec), np.max(delta_dec))

    if debug:
        key=raw_input("Enter any key to continue: ")

    xdata = delta_ra
    ydata = delta_dec

    print(xrange)
    print(yrange)

    # WISE W1 FWHM = 6.0 arcsecs
    radius = 6.0/2.0

    # include the psf radius so that sources that near edge are plotted
    itest = (xdata > (xrange[0] - radius)) & \
            (xdata < (xrange[1] + radius)) & \
            (ydata > (yrange[0] - radius)) & \
            (ydata < (yrange[1] + radius))

    xdata = xdata[itest]
    ydata = ydata[itest]
    ndata = len(xdata)

    if not overplot:
        plt.figure(figsize=(8,8))
        plt.axes().set_aspect('equal')

    patches = []
    ax = plt.gcf().gca()
    for x, y in zip(xdata, ydata):
        print(x, y, radius)
        circle = plt.Circle((x, y), radius, color='red', alpha=0.05)
        ax.add_artist(circle)

    plt.plot(xdata, ydata, '+', color='red')
    plt.xlim(xrange)
    plt.ylim(yrange)

    if showplot:
        plt.show()

    return plt
