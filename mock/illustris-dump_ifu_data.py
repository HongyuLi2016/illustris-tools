#!/usr/bin/env python
import numpy as np
import pyfits
from optparse import OptionParser
import sys


if __name__ == '__main__':
    parser = OptionParser()
    (options, args) = parser.parse_args()
    if len(args) != 1:
        print 'Error - please provide a folder name'
        sys.exit(1)
    path = args[0]

    spax_bin_data = np.loadtxt('{}/ifu/spaxel_bins.dat'.format(path))

    spax_bin_X = spax_bin_data[:, 0]
    spax_bin_Y = spax_bin_data[:, 1]
    spax_bin_id = spax_bin_data[:, 2].astype(int)

    bin_data = np.loadtxt('{}/ifu/voronoi_bins.dat'.format(path))

    bin_id = bin_data[:, 0].astype(int)
    bin_X = bin_data[:, 1]
    bin_Y = bin_data[:, 2]
    bin_area = bin_data[:, 3]
    bin_inuse = bin_data[:, 4].astype(int)

    ifu_data = np.loadtxt('{}/ifu/IFU_data'.format(path))
    bin_index = ifu_data[:, 0].astype(int)
    gh_v0 = ifu_data[:, 1]
    gh_v0_err = ifu_data[:, 2]
    gh_vd = ifu_data[:, 3]
    gh_vd_err = ifu_data[:, 4]
    gh_h3 = ifu_data[:, 5]
    gh_h3_err = ifu_data[:, 6]
    gh_h4 = ifu_data[:, 7]
    gh_h4_err = ifu_data[:, 8]
    v0 = ifu_data[:, 9]
    v0_err = ifu_data[:, 10]
    vd = ifu_data[:, 11]
    vd_err = ifu_data[:, 12]
    metal = ifu_data[:, 13]

    c1 = pyfits.Column(name='xbin', format='D', array=bin_X)
    c2 = pyfits.Column(name='ybin', format='D', array=bin_Y)
    c3 = pyfits.Column(name='v0', format='D', array=gh_v0)
    c4 = pyfits.Column(name='v0_err', format='D', array=gh_v0_err)
    c5 = pyfits.Column(name='vd', format='D', array=gh_vd)
    c6 = pyfits.Column(name='vd_err', format='D', array=gh_vd_err)
    c7 = pyfits.Column(name='h3', format='D', array=gh_h3)
    c8 = pyfits.Column(name='h3_err', format='D', array=gh_h3_err)
    c9 = pyfits.Column(name='h4', format='D', array=gh_h4)
    c10 = pyfits.Column(name='h4_err', format='D', array=gh_h4_err)
    c11 = pyfits.Column(name='metal', format='D', array=metal)

    c12 = pyfits.Column(name='rebin_x', format='D', array=spax_bin_X)
    c13 = pyfits.Column(name='rebin_y', format='D', array=spax_bin_Y)
    c14 = pyfits.Column(name='binid', format='D', array=spax_bin_id)

    coldefs1 = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11])
    coldefs2 = pyfits.ColDefs([c12, c13, c14])
    hdu = pyfits.PrimaryHDU()
    tbhdu1 = pyfits.BinTableHDU.from_columns(coldefs1)
    tbhdu2 = pyfits.BinTableHDU.from_columns(coldefs2)
    hdulist = pyfits.HDUList([hdu, tbhdu1, tbhdu2])
    hdulist.writeto('{}/ifu/IFU_gh.fits'.format(path), clobber=True)

    c1 = pyfits.Column(name='xbin', format='D', array=bin_X)
    c2 = pyfits.Column(name='ybin', format='D', array=bin_Y)
    c3 = pyfits.Column(name='v0', format='D', array=gh_v0)
    c4 = pyfits.Column(name='v0_err', format='D', array=gh_v0_err)
    c5 = pyfits.Column(name='vd', format='D', array=gh_vd)
    c6 = pyfits.Column(name='vd_err', format='D', array=gh_vd_err)
    c7 = pyfits.Column(name='metal', format='D', array=metal)

    c12 = pyfits.Column(name='rebin_x', format='D', array=spax_bin_X)
    c13 = pyfits.Column(name='rebin_y', format='D', array=spax_bin_Y)
    c14 = pyfits.Column(name='binid', format='D', array=spax_bin_id)

    coldefs1 = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7])
    coldefs2 = pyfits.ColDefs([c12, c13, c14])
    hdu = pyfits.PrimaryHDU()
    tbhdu1 = pyfits.BinTableHDU.from_columns(coldefs1)
    tbhdu2 = pyfits.BinTableHDU.from_columns(coldefs2)
    hdulist = pyfits.HDUList([hdu, tbhdu1, tbhdu2])
    hdulist.writeto('{}/ifu/IFU.fits'.format(path), clobber=True)

    # plt.plot(bin_X,bin_Y,'.r')
    # plt.plot(spax_bin_X,spax_bin_Y,'.k',alpha=0.1)
    # plt.show()
