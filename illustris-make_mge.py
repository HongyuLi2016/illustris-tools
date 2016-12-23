#!/usr/bin/env python
import numpy as np
# import matplotlib.pyplot as plt
# import pyfits
from optparse import OptionParser
# from glob import glob
import sys
from sectors_photometry import sectors_photometry
from mge_print_contours import mge_print_contours
from find_galaxy import find_galaxy
from mge_fit_sectors import mge_fit_sectors
import os
import util_illustris as ui


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-i', action='store_false', dest='iteration',
                      default=True, help='Iteration fit')
    (options, args) = parser.parse_args()
    if len(args) != 1:
        print 'Error - please provide a folder name'
        sys.exit(1)
    path = args[0]
    os.system('mkdir -p {}/mge'.format(path))
    try:
        img = np.load('{}/imgs/img_M.npy'.format(path))
    except:
        print 'Error - no image file found'
        sys.exit(1)
    # print img.shape
    part_numb = img/1.5e6*1e10
    ii = part_numb > 10
    # plt.imshow(np.log10(img))
    # img[~ii]=0.
    # plt.imshow(np.log10(img))
    # plt.show()
    pix2kpc = ui.scale_img   # 1 pixel = 0.25 kpc
    kpc2arcsec = ui.kpc2arcsec   # 1 kpc = 1.612 arcsec
    pix2arcsec = pix2kpc * kpc2arcsec   # pixel / arcsec
    level = img[ii].min()
    level *= 0.2
    old_img = img.copy()
    # img[55:85,115:145]=0
    lhy_f = find_galaxy(img, plot=1, fraction=0.1,  # level=level,
                        quiet=1, path='%s/mge/'%path)
    eps = lhy_f.eps
    theta = lhy_f.theta
    pa = np.mod(270 - theta, 180)
    xmed = lhy_f.xmed
    ymed = lhy_f.ymed
    xpeak = lhy_f.xpeak
    ypeak = lhy_f.ypeak
    lhy_s = sectors_photometry(img, eps, theta, xpeak, ypeak,
                               minlevel=level, plot=1, path='%s/mge/'%path)
    radius = lhy_s.radius
    angle = lhy_s.angle
    counts = lhy_s.counts
    qbound = 0.06
    lhy_mge = mge_fit_sectors(radius, angle, counts, eps, ngauss=15,
                              sigmaPSF=0, scale=pix2arcsec,
                              qbounds=[qbound, 0.999], linear=False, quiet=True,
                              outer_slope=4, bulge_disk=False,
                              plot=0, debug=False)
    sol = lhy_mge.sol.T
    absdev = lhy_mge.absdev

    if options.iteration:
        sol_old = sol
        absdev_old = absdev
        absdev_init = absdev
        while True:
            qbound += 0.05
            lhy_mge = mge_fit_sectors(radius, angle, counts, eps, ngauss=15,
                                      sigmaPSF=0, scale=pix2arcsec,
                                      qbounds=[qbound, 0.999], linear=False,
                                      quiet=True, outer_slope=4,
                                      bulge_disk=False, plot=0, debug=False)
            sol = lhy_mge.sol.T
            absdev = lhy_mge.absdev

            if (absdev/absdev_old > 1.10 or absdev/absdev_init >
                2.0 or qbound > 0.8):
                break
            absdev_old = absdev
            sol_old = sol
        sol = sol_old

    mge_print_contours(old_img, theta, xpeak, ypeak, sol.T,
                       sigmapsf=0, scale=pix2arcsec, magrange=5,
                       path='%s/mge/'%path)

    print 'total mass: %.4e'%(sol[:, 0].sum()*1e10)
    sol[:, 1] *= pix2arcsec
    sol[:, 0] = sol[:, 0]/(2*np.pi*sol[:, 1]**2*sol[:, 2])
    ml = 5.0  # mock stellar mass-to-light ratio
    # sol[:,0]*=1e10/(500.0**2)/ml
    sol[:, 0] *= 1e10/((pix2kpc*1e3)**2)/ml

    ff = open('%s/mge/mge.dat'%path, 'w+')
    print >>ff, 'Pa: %5.2f'%pa
    print >>ff, 'Eps: %5.3f'%eps
    print >>ff, 'Absdev: %4.3f'%absdev
    print >>ff, 'Xc Yc: %.3f %.3f'%(xmed, ymed)
    for i in range(len(sol[:, 0])):
        print >>ff, '%10.4e %10.2f %10.3f'%(sol[i, 0], sol[i, 1], sol[i, 2])
    ff.close()
    np.save('%s/mge/mge.npy'%path, [sol, pa, eps, xmed, ymed])
