#!/usr/bin/env python
import numpy as np
from optparse import OptionParser
import util_illustris as ui
import os
import sys
import pickle

parser = OptionParser()
(options, args) = parser.parse_args()
if len(args) != 1:
        print 'Error - please provide a folder name'
        sys.exit(1)

flist1 = np.genfromtxt('{}/merger/progenitor1/select.list'
                       .format(args[0]), dtype='S30')
flist2 = np.genfromtxt('{}/merger/progenitor2/select.list'
                       .format(args[0]), dtype='S30')
nsnap = flist1.shape[0]
mstar1 = np.zeros(nsnap)
mstar2 = np.zeros(nsnap)
mdark1 = np.zeros(nsnap)
mdark2 = np.zeros(nsnap)
pos1 = np.zeros([nsnap, 3])
pos2 = np.zeros([nsnap, 3])
vel1 = np.zeros([nsnap, 3])
vel2 = np.zeros([nsnap, 3])
orient_star1 = np.zeros([nsnap, 30, 3, 3])
orient_star2 = np.zeros([nsnap, 30, 3, 3])
orient_dark1 = np.zeros([nsnap, 100, 3, 3])
orient_dark2 = np.zeros([nsnap, 100, 3, 3])
a_star1 = np.zeros([nsnap, 3])
a_star2 = np.zeros([nsnap, 3])
c_star1 = np.zeros([nsnap, 3])
c_star2 = np.zeros([nsnap, 3])
a_dark1 = np.zeros([nsnap, 3])
a_dark2 = np.zeros([nsnap, 3])
c_dark1 = np.zeros([nsnap, 3])
c_dark2 = np.zeros([nsnap, 3])

for i in range(nsnap):
    z = ui.snap2z(int(flist1[i, 0][4:]))
    a = 1.0 / (1+z)
    with open('{}/merger/progenitor1/{}-{}/info.dat'
              .format(args[0], flist1[i, 0], flist1[i, 1])) as f:
        data1 = pickle.load(f)
    with open('{}/merger/progenitor2/{}-{}/info.dat'
              .format(args[0], flist2[i, 0], flist2[i, 1])) as f:
        data2 = pickle.load(f)
    with open('{}/merger/progenitor1/{}-{}/shape/shape.dat'
              .format(args[0], flist1[i, 0], flist1[i, 1])) as f:
        shape1 = pickle.load(f)
    with open('{}/merger/progenitor2/{}-{}/shape/shape.dat'
              .format(args[0], flist2[i, 0], flist2[i, 1])) as f:
        shape2 = pickle.load(f)
    mstar1[i] = data1['Subhalo']['SubhaloMassType'][4] * ui.massUnit
    mdark1[i] = data1['Subhalo']['SubhaloMassType'][1] * ui.massUnit
    pos1[i, :] = np.array(data1['Subhalo']['SubhaloPos']) * a / ui.h0
    vel1[i, :] = data1['Subhalo']['SubhaloVel']
    orient_star1[i, :, :, :] = shape1['eigenVectorsStar']
    orient_dark1[i, :, :, :] = shape1['eigenVectorsDark']

    mstar2[i] = data2['Subhalo']['SubhaloMassType'][4] * ui.massUnit
    mdark2[i] = data2['Subhalo']['SubhaloMassType'][1] * ui.massUnit
    pos2[i, :] = np.array(data2['Subhalo']['SubhaloPos']) * a / ui.h0
    vel2[i, :] = data2['Subhalo']['SubhaloVel']
    orient_star2[i, :, :, :] = shape2['eigenVectorsStar']
    orient_dark2[i, :, :, :] = shape2['eigenVectorsDark']
mratio_star = mstar1 / mstar2
mratio_dark = mdark1 / mdark2
relative_pos = pos1 - pos2
relative_vel = vel1 - vel2
L = np.cross(relative_pos, relative_vel, axis=1)
Lnorm = np.linalg.norm(L, axis=1)
for i in range(nsnap):
    L[i, :] /= Lnorm[i]
    Lnorm[i] *= (mstar1[i]*mstar2[i])/(mstar1[i]+mstar2[i])/1e10
    a_star1[i, :] = orient_star1[i, 11, 0, :]
    a_star2[i, :] = orient_star2[i, 11, 0, :]
    c_star1[i, :] = orient_star1[i, 11, 2, :]
    c_star2[i, :] = orient_star2[i, 11, 2, :]
    a_dark1[i, :] = orient_dark1[i, 39, 0, :]
    a_dark2[i, :] = orient_dark2[i, 39, 0, :]
    c_dark1[i, :] = orient_dark1[i, 39, 2, :]
    c_dark2[i, :] = orient_dark2[i, 39, 2, :]
rst = {}
rst['mratio_star'] = mratio_star
rst['mratio_dark'] = mratio_dark
rst['relative_pos'] = relative_pos
rst['relative_vel'] = relative_vel
rst['L'] = L
rst['Lnorm'] = Lnorm
rst['a_star1'] = a_star1
rst['a_star2'] = a_star2
rst['c_star1'] = c_star1
rst['c_star2'] = c_star2
rst['a_dark1'] = a_dark1
rst['a_dark2'] = a_dark2
rst['c_dark1'] = c_dark1
rst['c_dark2'] = c_dark2
with open('{}/merger/orbit.dat'.format(args[0]), 'wb') as f:
    pickle.dump(rst, f)
with open('{}/merger/orbit.txt'.format(args[0]), 'w') as f:
    print >> f, ('mratio_star:', mratio_star)
    print >> f, ('mratio_dark:', mratio_dark)
    print >> f, ('Lnorm', Lnorm)
    print >> f, ('L', L)
