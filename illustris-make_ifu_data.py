#!/usr/bin/env python
import os, sys, glob, math, time
from optparse import OptionParser

import numpy as np
# import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.spatial import cKDTree
from scipy.optimize import curve_fit
import util_illustris as ui


def inhullminmax(x, y, hull, minl, maxl, xmin, xmax, ymin, ymax):
    """Returns true if the point (x,y) is inside the convex hull
    and false otherwise"""

    R = math.sqrt(x * x + y * y)
    if R <= minl:
        return True
    if (x < xmin) or (x > xmax) or (y < ymin) or (y > ymax):
        return False

    for i in range(len(hull) - 1):
        x1 = hull[i+1, 0]
        y1 = hull[i+1, 1]
        x0 = hull[i, 0]
        y0 = hull[i, 1]

        dx = x1 - x0
        dy = y1 - y0

        if (y - y0) * dx - (x - x0) * dy > 0:
            return False

    return True


def hermite(n, x):
    """Calculate the value of the Hermite polynomial of degree n"""

    if n == 0:
        return 1.0
    elif n == 1:
        return x * math.sqrt(2.0)
    else:
        nd = float(n)
        term1 = math.sqrt(2.0 / nd) * x * hermite(n - 1, x)
        term2 = math.sqrt((nd - 1) / nd) * hermite(n - 2, x)
        return term1 - term2


def gh_series(v, meanv, sigma, h3, h4):
    """returns the value of a Gauss-Hermite series, used in estimating mean
    v, sigma, h3 and h4"""

    losvd = np.zeros(len(v))

    for i in xrange(len(v)):
        vnorm = (v[i] - meanv) / sigma

        Lvi = 1.0 + h3 * hermite(3, vnorm) + h4 * hermite(4, vnorm)

        losvd[i] = Lvi * math.exp(-0.5 * vnorm * vnorm) /\
            (sigma * math.sqrt(2.0 * math.pi))

    return losvd


def gh_gaussian(v, meanv, sigma):
    """returns the value of a Gaussian, used in estimating mean v, sigma"""

    losvd = np.zeros(len(v))

    for i in xrange(len(v)):
        vnorm = (v[i] - meanv) / sigma

        losvd[i] = math.exp(-0.5 * vnorm * vnorm) /\
            (sigma * math.sqrt(2.0 * math.pi))

    return losvd

if __name__ == '__main__':
    parser = OptionParser()
    (options, args) = parser.parse_args()
    if len(args) != 1:
        print 'Error - please provide a folder name'
        sys.exit(1)
    path = args[0]
    os.system('mkdir -p %s/ifu/kenimatic_fit'%path)
    # read hull_file
    with open('%s/ifu/IFU_hull'%path, 'r') as f:
        ctrl_fields = f.readline().split()
    num_hull_points = int(ctrl_fields[0])
    hull_minR = float(ctrl_fields[1])
    hull_maxR = float(ctrl_fields[2])
    hull_xmin = float(ctrl_fields[3])
    hull_xmax = float(ctrl_fields[4])
    hull_ymin = float(ctrl_fields[5])
    hull_ymax = float(ctrl_fields[6])

    print
    print 'Number of hull points = {0}'.format(num_hull_points)
    hull_points = np.loadtxt('%s/ifu/IFU_hull'%path, skiprows=1)
    hullx = hull_points[:, 0]
    hully = hull_points[:, 1]

    # read ifu bin file
    bin_data = np.loadtxt('%s/ifu/voronoi_bins.dat'%path)

    bin_id    = bin_data[:, 0].astype(int)
    bin_X     = bin_data[:, 1]
    bin_Y     = bin_data[:, 2]
    bin_area  = bin_data[:, 3]
    bin_inuse = bin_data[:, 4].astype(int)

    num_ifu_bins = len(bin_id)
    print 'Number of IFU bins = {0}'.format(num_ifu_bins)

    XYvalues = np.zeros((num_ifu_bins, 2))
    XYvalues[:, 0] = bin_X[:]
    XYvalues[:, 1] = bin_Y[:]

    # create a kd tree for binning particles
    ifu_tree = cKDTree(XYvalues)

    # read particle coordinates
    data = np.load('%s/imgs/coordinates_star.npy'%path)
    X = data[:, 0] * ui.kpc2arcsec  # convert position unit to arcsec
    Y = data[:, 1] * ui.kpc2arcsec
    Z = data[:, 2] * ui.kpc2arcsec
    VX = data[:, 3]
    VY = data[:, 4]
    VZ = data[:, 5]
    M  = data[:, 6]
    L  = data[:, 7]
    Metal = data[:, 8]
    # what = np.histogram2d(X,Z,range=[[-20,20],[-20,20]],bins=200)
    # img = what[0]
    # print img.shape
    # plt.imshow(np.log10(img.T))
    # plt.show()
    num_particles = len(X)
    print 'Number of star particles = {0}'.format(num_particles)

    particle_bins = np.zeros(num_particles, dtype=int)
    for i in xrange(num_particles):
        # bin assuming los is Y axis
        if inhullminmax(X[i], Z[i], hull_points, hull_minR, hull_maxR,
                        hull_xmin, hull_xmax, hull_ymin, hull_ymax):
            data_point = [X[i], Z[i]]
            distance, bindex = ifu_tree.query(data_point)
            particle_bins[i] = bindex
        else:
            particle_bins[i] = -1

    inbin_mask = particle_bins >= 0

    inbin_bin  = particle_bins[inbin_mask]
    inbin_VY   = VY[inbin_mask]
    inbin_M    = M[inbin_mask]
    inbin_L    = L[inbin_mask]
    inbin_Metal= Metal[inbin_mask]
    print 'Number of binned particles = {0}'.format(len(inbin_bin))

    IFU_data = open('%s/ifu/IFU_data'%path, 'w')

    for bin_index in xrange(num_ifu_bins):
        this_bin_mask = inbin_bin == bin_index

        # create the velocity histogram

        vrange = [-800.0, 800.0]
        inbin_velocity = inbin_VY[this_bin_mask]
        inbin_Mwt      = inbin_M[this_bin_mask]
        hist, bin_edges = np.histogram(inbin_velocity, range=vrange, bins=50,
                                       normed=True, weights=inbin_Mwt)
        edges = np.array([0.5 * (bin_edges[i+1] + bin_edges[i])
                          for i in xrange(len(bin_edges) - 1)])
        edgesplot = edges.copy()
        # remove zero points

        zero_mask = hist != 0.0
        hist  = hist[zero_mask]
        edges = edges[zero_mask]

        # fit G-H series to the histogrammed data using curve_fit from scipy
        # LHY variables with suffix g is add by lhy to fit a pure gaussian LOSVD

        p0 = [np.mean(edges), np.std(edges), 0.0, 0.0]
        p0g = [np.mean(edges), np.std(edges)]
        try:
            popt, pcov = curve_fit(gh_series, edges, hist, p0, maxfev=1500)
            poptg, pcovg = curve_fit(gh_gaussian, edges, hist, p0g, maxfev=1500)
        except RuntimeError:
            print 'Error - failed processing bin {0}, continuing but check output'.format(bin_index)
            continue
        else:
            perr = np.sqrt(np.diag(pcov))
            perrg = np.sqrt(np.diag(pcovg))

        Lv = gh_series(edgesplot, popt[0], popt[1], popt[2], popt[3])
        Lvg = gh_gaussian(edgesplot, poptg[0], poptg[1])

        plt.title('Subhalo = {0} - los velocity distribution, bin = {1}'
                  .format(path, bin_index))
        plt.xlabel('Los velocity (km/s)')
        plt.plot(edges, hist, 'bo', edgesplot, Lv, 'r-', edgesplot, Lvg, 'c-')
        ymin, ymax = plt.ylim()
        plt.vlines(popt[0], 0.0, ymax)
        sigma_lines = [popt[0] - popt[1], popt[0] + popt[1]]
        plt.vlines(sigma_lines, 0.0, ymax, linestyles='dashed')
        plt.xlim(vrange)
        plt.savefig('{0}/ifu/kenimatic_fit/zzz_{1:03d}'.format(path, bin_index))
        plt.clf()

        # calculate metalicity
        inthisbin_Metal  = inbin_Metal[this_bin_mask]
        massMetal = (inthisbin_Metal * inbin_Mwt).sum()
        IFU_data.write('{:3d} {:+e} {:+e} {:+e} {:+e} {:+e} {:+e} {:+e}'
                       ' {:+e} {:+e} {:+e} {:+e} {:+e} {:+e}\n'.format(
                           bin_index, popt[0], perr[0], popt[1], perr[1],
                           popt[2], perr[2], popt[3], perr[3], poptg[0],
                           perrg[0], poptg[1], perrg[1], massMetal))

    IFU_data.close()
