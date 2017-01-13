#!/usr/bin/env python
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
import pyfits
from optparse import OptionParser
from JAM.utils.velocity_plot import velocity_plot
from scipy import stats


def v_plot(path):
    if not os.path.isfile('%s/ifu/IFU.fits'%path):
        print 'No IFU.fits in %s'%path
        exit(1)
    hdulist = pyfits.open('%s/ifu/IFU.fits'%path)[1]
    tem = hdulist.data
    x0 = tem['xbin']
    y0 = tem['ybin']
    v0 = tem['v0']
    v0_err = tem['v0_err']
    vd = tem['vd']
    vd_err = tem['vd_err']
    r = np.sqrt(x0**2+y0**2)
    ii = np.where(r < 3.0)
    mv0 = np.mean(v0[ii])
    vel = v0-mv0
    good = (vel**2+vd**2)**0.5 < 800.
    x0 = x0[good]
    y0 = y0[good]
    vel = vel[good]
    v0_err = v0_err[good]
    vd = vd[good]
    vd_err = vd_err[good]

    fig = plt.figure(figsize=(4*1.5, 3.3*1.5))
    fig.subplots_adjust(left=0.07, bottom=0.05, right=0.88,
                        top=0.98, wspace=0.6, hspace=0.01)

    ax1 = fig.add_subplot(2, 2, 1)
    nans = np.isnan(vel)
    vmax = stats.scoreatpercentile(vel[~nans], 98.0)
    norm = colors.Normalize(vmin=-vmax, vmax=vmax)
    velocity_plot(x0, y0, vel, markersize=0.2, norm=norm,
                  ax=ax1, text='$\mathbf{V^{*}}$', equal=True,
                  xreverse=False)
    ax2 = fig.add_subplot(2, 2, 2)
    velocity_plot(x0, y0, v0_err.clip(0., 100.), markersize=0.2,
                  ax=ax2, text='$\mathbf{V^{*}_{err}}$', equal=True,
                  xreverse=False)
    ax3 = fig.add_subplot(2, 2, 3)
    velocity_plot(x0, y0, vd.clip(0, 600.), markersize=0.2,
                  ax=ax3, text='$\mathbf{\sigma^{*}}$',
                  xreverse=False)
    ax4 = fig.add_subplot(2, 2, 4)
    velocity_plot(x0, y0, vd_err.clip(0, 100.), markersize=0.2,
                  ax=ax4, text='$\mathbf{\sigma^{*}_{err}}$',
                  xreverse=False)
    fig.savefig('%s/ifu/IFU.png'%path, dpi=300)

if __name__ == '__main__':
    parser = OptionParser()
    (options, args) = parser.parse_args()
    if len(args) != 1:
        print 'Error - please provide a folder name'
        sys.exit(1)
    path = args[0]
    v_plot(path)
