#!/usr/bin/env python
import numpy as np
import util_illustris as ui
from optparse import OptionParser
import sys
import os
import pickle


def shape(R, xpart):
    mpart = np.ones(xpart.shape[0])
    axisRatios = np.zeros([len(R), 2])
    eigenVectors = np.zeros([len(R), 3, 3])
    for i in range(len(R)):
        ba, ca, angle, Tiv = ui.get_shape(xpart, mpart, Rb=R[i])
        for j in range(3):
            if Tiv[j, 2] < 0:
                Tiv[j, :] *= -1
        # eigen vectors Tiv[0, :] Tiv[1, :] Tiv[2, :]
        axisRatios[i, :] = np.array([ba, ca])
        eigenVectors[i, :, :] = Tiv
    return axisRatios, eigenVectors


def get_shape_main(path):
    os.system('mkdir -p {}/shape'.format(path))
    with open('{}/info.dat'.format(path), 'rb') as f:
        info = pickle.load(f)
    snapNum = int(info['SnapshotNumber'])
    z = ui.snap2z(snapNum)
    pos = info['Subhalo']['SubhaloPos'] / (1 + z) / ui.h0

    # half mass radius [kpc] (physical unit)
    hmr_star = info['Subhalo']['SubhaloHalfmassRadType'][4] / (1 + z) / ui.h0
    hmr_dark = info['Subhalo']['SubhaloHalfmassRadType'][1] / (1 + z) / ui.h0
    xpart_star = ui.read_cutout('{}/cutout.hdf5'.format(path), PartType=4,
                                key='Coordinates', z=z) - pos
    xpart_dark = ui.read_cutout('{}/cutout.hdf5'.format(path), PartType=1,
                                key='Coordinates', z=z) - pos
    Rstar = np.linspace(3.0, 2.5*hmr_star, 30)
    Rdark = np.linspace(3.0, 2.5*hmr_dark, 100)
    axisRatiosStar, eigenVectorsStar = shape(Rstar, xpart_star)
    axisRatiosDark, eigenVectorsDark = shape(Rdark, xpart_dark)
    rst = {'Rstar': Rstar, 'hmr_star': hmr_star, 'axisRatiosStar':
           axisRatiosStar, 'eigenVectorsStar': eigenVectorsStar,
           'Rdark': Rdark, 'hmr_dark': hmr_dark, 'axisRatiosDark':
           axisRatiosDark, 'eigenVectorsDark': eigenVectorsDark}
    with open('{}/shape/shape.dat'.format(path), 'wb') as f:
        pickle.dump(rst, f)


if __name__ == '__main__':
    parser = OptionParser()
    (options, args) = parser.parse_args()
    if len(args) != 1:
        print 'Error - please provide a folder name'
        sys.exit(1)
    get_shape_main(args[0])
