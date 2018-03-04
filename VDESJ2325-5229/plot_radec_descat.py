import sys

sys.path.append('/home/rgm/soft/python/lib/')
from librgm.plotid import plotid


def plot_radec_descat(data=None,
                      source=None,
                      radec_centre=None,
                      xrange = [-2.5, 2.5],
                      yrange = [-2.5, 2.5],
                      wavebands=['g', 'r', 'i', 'z', 'y'],
                      coadd=True, multiBand=True,
                      singleEpoch=False,
                      showplot=False,
                      debug=True):
    """


    """

    import numpy as np

    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Ellipse

    colors=['blue','green','orange','red', 'maroon']

    print('filename:', data.meta['filename'])

    print('Number of rows:', len(data))
    if debug:
        data.info()
        data.info('stats')
        print('zoom:', zoom)
        print('filename:', data.meta['filename'])
        print('ra_centre', radec_centre[0])
        print('dec_centre', radec_centre[1])
        print('xrange:', xrange)
        print('yrange:', yrange)

    itest = np.unique(data['TILENAME'])
    print('Tiles:', itest)

    ra = data['RA']
    dec = data['DEC']

    COADD_OBJECTS_ID = data['COADD_OBJECTS_ID']

    ra_min = np.min(ra)
    ra_max = np.max(ra)
    # convert to arc seconds
    delta_ra = (ra - radec_centre[0]) * 3600.0 * \
               np.cos(np.deg2rad(radec_centre[1]))

    if debug:
        print('ra_min:', ra_min)
        print('ra_max:', ra_max)
        print('Delta RA range:',
              np.min(delta_ra), np.max(delta_ra))

    dec_min = np.min(dec)
    dec_max = np.max(dec)
    delta_dec = (dec - radec_centre[1]) * 3600.0

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
    itest = (xdata > xrange[0]) & (xdata < xrange[1]) & \
            (ydata > yrange[0]) & (ydata < yrange[1])

    xdata = xdata[itest]
    ydata = ydata[itest]
    ndata = len(xdata)

    plt.figure(figsize=(8,8))
    plt.axes().set_aspect('equal')
    plt.plot(xdata, ydata, '.k', label='COADD: '+str(ndata))

    plt.xlim(xrange)
    plt.ylim(yrange)

    infile = data.meta['filename']
    plt.title(infile, fontsize='medium')
    plt.grid(True)

    WAVEBANDS = ['G','R','I','Z','Y']

    print("id, object_id, WAVEBAND, NEPOCHS, MAG_AUTO, " + \
              "Dra, Ddec, " + \
              "FLUX_RADIUS, KRON_RADIUS, " + \
              "A, B, PA, Aspect Ratio(A/B)")

    iband = -1
    for WAVEBAND in WAVEBANDS:
        iband = iband + 1
        ra = data['ALPHAWIN_J2000_'+WAVEBAND]
        dec = data['DELTAWIN_J2000_'+WAVEBAND]

        # used to convert pixels to arc seconds
        PIXEL_SIZE=0.27

        delta_ra = (ra - radec_centre[0])*3600.0 * \
                   np.cos(np.deg2rad(radec_centre[1]))
        delta_dec = (dec - radec_centre[1])*3600.0

        xdata= delta_ra
        ydata= delta_dec

        # limit the data
        itest = (xdata > xrange[0]) & (xdata < xrange[1]) & \
            (ydata > yrange[0]) & (ydata < yrange[1])

        xdata = xdata[itest]
        ydata = ydata[itest]
        ndata = len(xdata)

        delta_ra = delta_ra[itest]
        delta_dec = delta_dec[itest]

        # maybe combine this with the assignment above
        COADD_OBJECTS_ID = data['COADD_OBJECTS_ID'][itest]
        OBJECT_NUMBER = data['OBJECT_NUMBER'][itest]
        RADEC_SOURCE_BAND = data['RADEC_SOURCE_BAND'][itest]

        # Detection image
        A_IMAGE = data['A_IMAGE'][itest]*PIXEL_SIZE
        B_IMAGE = data['B_IMAGE'][itest]*PIXEL_SIZE
        THETA_IMAGE = data['THETA_IMAGE'][itest]

        KRON_RADIUS = data['KRON_RADIUS'][itest]*PIXEL_SIZE

        XMIN_IMAGE=data['XMIN_IMAGE'][itest]
        XMAX_IMAGE=data['XMAX_IMAGE'][itest]
        YMIN_IMAGE=data['YMIN_IMAGE'][itest]
        YMAX_IMAGE=data['YMAX_IMAGE'][itest]

        # Measurement image
        ALPHAWIN_J2000 = data['ALPHAWIN_J2000_'+WAVEBAND][itest]
        DELTAWIN_J2000 = data['DELTAWIN_J2000_'+WAVEBAND][itest]
        XWIN_IMAGE = data['XWIN_IMAGE_'+WAVEBAND][itest]
        YWIN_IMAGE = data['YWIN_IMAGE_'+WAVEBAND][itest]

        NEPOCHS = data['NEPOCHS_'+WAVEBAND][itest]
        MAG_AUTO = data['MAG_AUTO_'+WAVEBAND][itest]
        FLUX_RADIUS = data['FLUX_RADIUS_'+WAVEBAND][itest]*PIXEL_SIZE

        AWIN_IMAGE = data['AWIN_IMAGE_'+WAVEBAND][itest]*PIXEL_SIZE
        BWIN_IMAGE = data['BWIN_IMAGE_'+WAVEBAND][itest]*PIXEL_SIZE
        THETAWIN_IMAGE = data['THETAWIN_IMAGE_'+WAVEBAND][itest]

        ELLIP1MODEL_WORLD = data['ELLIP1MODEL_WORLD_'+WAVEBAND][itest]
        ELLIP2MODEL_WORLD = data['ELLIP2MODEL_WORLD_'+WAVEBAND][itest]

        FWHM_WORLD = data['FWHM_WORLD_'+WAVEBAND][itest]

        ISOAREA_WORLD = data['ISOAREA_WORLD_'+WAVEBAND][itest]

        FLAGS = data['FLAGS_'+WAVEBAND][itest]

        alpha=0.1
        i =-1

        for id in xdata:
            i = i + 1

            circle = Circle([delta_ra[i], delta_dec[i]], 0.25,
                edgecolor='none', facecolor=colors[iband], alpha=alpha)
            plt.gca().add_patch(circle)

            width = AWIN_IMAGE[i]
            height = BWIN_IMAGE[i]
            angle = THETAWIN_IMAGE[i] + 90.0
            coadd_objects_id = COADD_OBJECTS_ID[i]

            print(i, coadd_objects_id,
                OBJECT_NUMBER[i],
                WAVEBAND,
                RADEC_SOURCE_BAND[i],
                "{:4d}".format(NEPOCHS[i]),
                "{:8.2f}".format(MAG_AUTO[i]),
                "{:6.2f}".format(delta_ra[i]),
                "{:6.2f}".format(delta_dec[i]),
                "{:6.2f}".format(FLUX_RADIUS[i]),
                "{:6.2f}".format(KRON_RADIUS[i]),
                "{:6.2f}".format(width),
                "{:6.2f}".format(height),
                "{:7.1f}".format(angle),
                "{:6.3f}".format(width/height))

            print(i, coadd_objects_id,
                OBJECT_NUMBER[i],
                WAVEBAND,
                "{:8.2f}".format(A_IMAGE[i]),
                "{:8.2f}".format(B_IMAGE[i]),
                "{:7.1f}".format(THETA_IMAGE[i]))

            print(i, coadd_objects_id,
                OBJECT_NUMBER[i],
                WAVEBAND,
                "{:7.1f}".format(XWIN_IMAGE[i]),
                "{:7.1f}".format(YWIN_IMAGE[i]),
                "{:7.1f}".format(XMIN_IMAGE[i]),
                "{:7.1f}".format(XMAX_IMAGE[i]),
                "{:7.1f}".format(YMIN_IMAGE[i]),
                "{:7.1f}".format(YMAX_IMAGE[i]),
                "{:5.1f}".format(XMAX_IMAGE[i]-XMIN_IMAGE[i]),
                "{:5.1f}".format(YMAX_IMAGE[i]-YMIN_IMAGE[i]),
                "{:7.1f}".format((XMIN_IMAGE[i]+XMAX_IMAGE[i])/2.0),
                "{:7.1f}".format((YMIN_IMAGE[i]+YMAX_IMAGE[i])/2.0),
                "{:4d}".format(FLAGS[i]))

            ellipse = Ellipse([delta_ra[i], delta_dec[i]],
                width=width/2.0, height=height/2.0, angle=angle,
                edgecolor='none', facecolor=colors[iband], alpha=alpha)
            plt.gca().add_patch(ellipse)

        plt.plot(delta_ra, delta_dec, '.', color=colors[iband],
            label=WAVEBAND + ': ' + str(ndata))

    plt.suptitle('Windowed positions')

    plt.xlabel('Delta RA (arc seconds)')
    plt.ylabel('Delta Dec (arc seconds)')
    plt.legend(fontsize='medium')
    plotid()

    plotfile = source + '_COADD_radec.png'
    plt.savefig(plotfile)
    #plt.clf()
    print('Saving: ', plotfile)

    if showplot:
        plt.show()

    return plt
