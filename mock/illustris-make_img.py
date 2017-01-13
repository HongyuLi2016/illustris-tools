#!/usr/bin/env python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# import pyfits
from optparse import OptionParser
from glob import glob
import os
import sys
import util_illustris as ui
import warnings

warnings.simplefilter("ignore")


def make_img(x, L, scale=0.5, boxsize=50.0):
    '''
    make a mock image, pixle size = 0.5 kpc/pix,
    image size = 50.0 kpc * 50.0 kpc
    '''
    ngrid = int(boxsize/scale)
    img = np.zeros([ngrid, ngrid])
    for n in range(len(x)):
        # assume y is LOS
        i = int(np.rint(x[n, 0]/scale) + int(ngrid/2))
        k = int(np.rint(x[n, 2]/scale) + int(ngrid/2))
        if i >= 0 and i < ngrid and k >= 0 and k < ngrid:
            # axis=1 -> x  axis=0 ->z
            img[k, i] += L[n]
    return img


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-i', action='store', type='float', dest='i',
                      default=45.0, help='inclination')
    parser.add_option('-p', action='store', type='float', dest='phi',
                      default=0.0, help='phi')
    (options, args) = parser.parse_args()
    if len(args) != 1:
        print 'Error - please provide a folder name'
        sys.exit(1)
    path = args[0]
    cutoutName = glob('{}/cutout*.hdf5'.format(path))
    if len(cutoutName) != 1:
        print 'Error - No hdf5 particle file in {}'.format(path)

    # read particle coordinates
    xpart_star = ui.read_cutout(cutoutName[0], PartType=4,
                                key='Coordinates')
    vpart_star = ui.read_cutout(cutoutName[0], PartType=4,
                                key='Velocities')
    mpart_star = ui.read_cutout(cutoutName[0], PartType=4,
                                key='Masses')
    Meta = ui.read_cutout(cutoutName[0], PartType=4,
                          key='GFM_Metallicity')
    lpart = ui.read_cutout(cutoutName[0], PartType=4,
                           key='GFM_StellarPhotometrics')
    lrpart = 10**((4.58 - lpart[:, 5])/2.5)/1e10  # r band luminosity
    xcenter_star = ui.find_center(xpart_star, mpart=mpart_star, percent=20.0)
    xpart_star -= xcenter_star
    vsys_star = ui.vel_center(xpart_star, vpart_star, mpart=mpart_star, R=10.0)
    vpart_star -= vsys_star
    '''
    deltaD = 10.0
    old_massCenterStar = np.average(xpart_old_star, axis=0, weights=mpart_star)
    xpart_old_star -= old_massCenterStar
    # plt.hist(xpart_old_star[:,0],bins=100)
    # plt.hist(xpart_old_star[:,1],bins=100)
    # plt.hist(xpart_old_star[:,2],bins=100)
    # plt.show()

    # iteratively find out the mass center
    while deltaD > 0.01:
        rstar = (xpart_old_star[:, 0]**2 + xpart_old_star[:, 1]**2 +
                 xpart_old_star[:, 2]**2)**0.5
        iInStar = rstar < np.percentile(rstar, 20.0)
        massCenterStar = np.average(xpart_old_star[iInStar, :],
                                    axis=0, weights=mpart_star[iInStar])
        # print massCenterStar,old_massCenterStar
        deltaD = np.sum((massCenterStar - old_massCenterStar)**2)**0.5
        # print deltaD
        old_massCenterStar = massCenterStar.copy()
        xpart_old_star -= massCenterStar
        # plt.hist(xpart_old_star[:,0],bins=100)
        # plt.hist(xpart_old_star[:,1],bins=100)
        # plt.hist(xpart_old_star[:,2],bins=100)
        # plt.show()
    '''

    xpart_dark = ui.read_cutout(cutoutName[0], PartType=1, key='Coordinates')
    vpart_dark = ui.read_cutout(cutoutName[0], PartType=1, key='Velocities')
    xcenter_dark = ui.find_center(xpart_dark, percent=20.0)
    xpart_dark -= xcenter_dark
    vsys_dark = ui.vel_center(xpart_dark, vpart_dark, R=10.0)
    vpart_dark -= vsys_dark

    # difference between halo and galaxy
    centerDiff = xcenter_dark - xcenter_star
    vsysDiff = vsys_dark - vsys_star

    # get galaxy shape and axis
    ba, ca, angle, Tiv = ui.get_shape(xpart_star, mpart_star, Rb=10.)

    # coordinates in 3-axis-system
    # positions and velocities in 3 axis coordinate system
    xpart_axis = np.dot(Tiv, xpart_star.T).T
    vpart_axis = np.dot(Tiv, vpart_star.T).T

    # rotate by the phi alone z axis and then by the inclination
    # alone new x axis
    # left hand system, with y point to LOS when edge on
    phi = options.phi / 180.0 * np.pi
    inc = (90.0 - options.i) / 180.0 * np.pi  # i =90 => edge-on
    M_rot_phi = np.array([[np.cos(phi), -np.sin(phi), 0],
                          [np.sin(phi), np.cos(phi), 0], [0, 0, 1]])
    M_rot_inc = np.array([[1, 0, 0], [0, np.cos(inc), -np.sin(inc)],
                          [0, np.sin(inc), np.cos(inc)]])
    M_rot = np.dot(M_rot_inc, M_rot_phi)
    # positions and velocities in mock obs system
    xpart = np.dot(M_rot, xpart_axis.T).T
    vpart = np.dot(M_rot, vpart_axis.T).T

    boxsize_img = ui.boxsize_img
    scale_img = ui.scale_img
    img_M = make_img(xpart, mpart_star, scale=scale_img, boxsize=boxsize_img)
    img_L = make_img(xpart, lrpart, scale=scale_img, boxsize=boxsize_img)
    boxsize_ifu = ui.boxsize_ifu
    scale_ifu = ui.scale_ifu
    img_ifu = make_img(xpart, lrpart*0.0+1, scale=scale_ifu,
                       boxsize=boxsize_ifu)
    # print lrpart.sum(),mpart_star.sum()
    # save necessary files
    os.system('mkdir -p {}/imgs'.format(path))

    np.save('{}/imgs/img_M.npy'.format(path), img_M)
    np.save('{}/imgs/img_L.npy'.format(path), img_L)
    np.save('{}/imgs/img_ifu.npy'.format(path), img_ifu)

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.imshow(np.log10(img_M), origin='lower',
              extent=(-boxsize_img/2, boxsize_img/2, -boxsize_img/2,
                      boxsize_img/2))
    ax.set_xlabel('kpc')
    ax.set_ylabel('kpc')
    fig.savefig('{}/imgs/img_M.png'.format(path), dpi=300)

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.imshow(np.log10(img_L), origin='lower',
              extent=(-boxsize_img/2, boxsize_img/2, -boxsize_img/2,
                      boxsize_img/2))
    ax.set_xlabel('kpc')
    ax.set_ylabel('kpc')
    fig.savefig('{}/imgs/img_L.png'.format(path), dpi=300)

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.imshow(np.log10(img_ifu), origin='lower',
              extent=(-boxsize_ifu/2, boxsize_ifu/2, -boxsize_ifu/2,
                      boxsize_ifu/2))
    ax.set_xlabel('kpc')
    ax.set_ylabel('kpc')
    fig.savefig('{}/imgs/img_ifu.png'.format(path), dpi=300)

    with open('{}/imgs/make_img.log'.format(path), 'w') as ff:
        print >>ff, 'inclination: {:.2f}'.format(options.i)
        print >>ff, 'phi: {:.2f}'.format(options.phi)
        print >>ff, 'boxsize_img: {:.2f} kpc'.format(boxsize_img)
        print >>ff, 'boxsize_ifu: {:.2f} kpc'.format(boxsize_ifu)
        print >>ff, 'scale_img: {:.2f} kpc/pixel'.format(scale_img)
        print >>ff, 'scale_ifu: {:.2f} kpc/pixel'.format(scale_ifu)
        print >>ff, 'Stellar mass: {:.4e} M_solar'.format(mpart_star.sum()*1e10)
        print >>ff, ('r band luminosity: {:.4e} L_solar'
                     .format(lrpart.sum()*1e10))
        print >>ff, ('Averaged M*/L: {:.3f} M_solar/L_solar'
                     .format(mpart_star.sum()/lrpart.sum()))
        print >>ff, ('Mass centre offset: {:.2f} {:.2f} {:.2f} kpc'
                     .format(centerDiff[0], centerDiff[1], centerDiff[2]))
        print >>ff, ('System velocity offset: {:.2f} {:.2f} {:.2f} km/s'
                     .format(vsysDiff[0], vsysDiff[1], vsysDiff[2]))
        print >>ff, 'ba={:.2f}  ca={:.2f}'.format(ba, ca)
    # x in kpc, v in km/s, mpart_star in 10^10 M_solar, lrpart in 10^10 L_solar
    data = np.zeros([xpart.shape[0], 9])
    data[:, 0:3] = xpart
    data[:, 3:6] = vpart
    data[:, 6] = mpart_star
    data[:, 7] = lrpart
    data[:, 8] = Meta
    np.save('{}/imgs/coordinates_star.npy'.format(path), data)
