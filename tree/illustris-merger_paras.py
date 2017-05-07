#!/usr/bin/env python
import numpy as np
from optparse import OptionParser
import util_illustris as ui
import os
import sys
import pickle

def get_angle(vect1, vect2, acute=True):
    '''
    Calculate the angle between 2 vectors
    '''
    vect1 = vect1 / np.linalg.norm(vect1)
    vect2 = vect2 / np.linalg.norm(vect2)
    cosi = np.dot(vect1, vect2)
    angle = np.degrees(np.arccos(cosi))
    if acute:
        if angle >= 90.0:
            angle = 180 - angle
    return angle

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
spinValue_star1 = np.zeros([nsnap, 30])
spinValue_star2 = np.zeros([nsnap, 30])
spinVector_star1 = np.zeros([nsnap, 30, 3])
spinVector_star2 = np.zeros([nsnap, 30, 3])
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
spinVector1 = np.zeros([nsnap, 3])
spinVector2 = np.zeros([nsnap, 3])
spinValue1 = np.zeros(nsnap)
spinValue2 = np.zeros(nsnap)
orbit_plane = np.zeros([nsnap-1, 3])
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
    with open('{}/merger/progenitor1/{}-{}/spin/spin.dat'
              .format(args[0], flist1[i, 0], flist1[i, 1])) as f:
        spin1 = pickle.load(f)
    with open('{}/merger/progenitor2/{}-{}/spin/spin.dat'
              .format(args[0], flist2[i, 0], flist2[i, 1])) as f:
        spin2 = pickle.load(f)
    cutout1 = ui.cutout('{}/merger/progenitor1/{}-{}/cutout.hdf5'
                        .format(args[0], flist1[i, 0], flist1[i, 1]))
    cutout2 = ui.cutout('{}/merger/progenitor2/{}-{}/cutout.hdf5'
                        .format(args[0], flist2[i, 0], flist2[i, 1]))
    mstar1[i] = data1['Subhalo']['SubhaloMassType'][4] * ui.massUnit
    mdark1[i] = data1['Subhalo']['SubhaloMassType'][1] * ui.massUnit
    pos1[i, :] = np.array(data1['Subhalo']['SubhaloPos'])
    hmr1 = data1['Subhalo']['SubhaloHalfmassRadType'][4]

    vpart_star1 = cutout1.vpart_star * np.sqrt(a)
    xpart_star1 = cutout1.xpart_star
    mpart_star1 = cutout1.mass_star
    vcenter_star1 = ui.vel_center(xpart_star1-pos1[i, :], vpart_star1,
                                  mpart=mpart_star1, R=hmr1*0.5)
    vel1[i, :] = vcenter_star1
    orient_star1[i, :, :, :] = shape1['eigenVectorsStar']
    orient_dark1[i, :, :, :] = shape1['eigenVectorsDark']
    spinValue_star1[i, :] = spin1['spinValueStar']
    spinVector_star1[i, :, :] = spin1['spinVectorsStar']

    mstar2[i] = data2['Subhalo']['SubhaloMassType'][4] * ui.massUnit
    mdark2[i] = data2['Subhalo']['SubhaloMassType'][1] * ui.massUnit
    pos2[i, :] = np.array(data2['Subhalo']['SubhaloPos'])
    hmr2 = data2['Subhalo']['SubhaloHalfmassRadType'][4]

    vpart_star2 = cutout2.vpart_star * np.sqrt(a)
    xpart_star2 = cutout2.xpart_star
    mpart_star2 = cutout2.mass_star
    vcenter_star2 = ui.vel_center(xpart_star2-pos2[i, :], vpart_star2,
                                  mpart=mpart_star2, R=hmr2*0.5)
    vel2[i, :] = vcenter_star2
    orient_star2[i, :, :, :] = shape2['eigenVectorsStar']
    orient_dark2[i, :, :, :] = shape2['eigenVectorsDark']
    spinValue_star2[i, :] = spin2['spinValueStar']
    spinVector_star2[i, :, :] = spin2['spinVectorsStar']
mratio_star = mstar1 / mstar2
mratio_dark = mdark1 / mdark2
relative_pos = pos1 - pos2
for i in range(relative_pos.shape[0]):
    for j in range(relative_pos.shape[1]):
        if relative_pos[i, j] > (75*1e3 * 0.6):
            relative_pos[i, j] -= 75 * 1e3
        if relative_pos[i, j] < (-75*1e3 * 0.6):
            relative_pos[i, j] += 75 * 1e3
relative_pos *= a / ui.h0
relative_vel = vel1 - vel2
for i in range(nsnap-1):
    tem = np.cross(relative_pos[i, :], relative_pos[i+1, :])
    orbit_plane[i, :] = tem / np.linalg.norm(tem)
r = np.linalg.norm(relative_pos, axis=1)
L = np.cross(relative_pos, relative_vel, axis=1)
Lnorm = np.linalg.norm(L, axis=1)
for i in range(nsnap):
    L[i, :] /= Lnorm[i]
    Lnorm[i] *= (mstar1[i]*mstar2[i])/(mstar1[i]+mstar2[i])/1e10
    spinValue1[i] = spinValue_star1[i, 11]
    spinVector1[i, :] = spinVector_star1[i, 11, :]
    spinValue2[i] = spinValue_star2[i, 11]
    spinVector2[i, :] = spinVector_star2[i, 11, :]
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
rst['orbit_plane'] = orbit_plane
rst['spinValue1'] = spinValue1
rst['spinVector1'] = spinVector1
rst['spinValue2'] = spinValue2
rst['spinVector2'] = spinVector2
cs12 = np.zeros(nsnap)
as12 = np.zeros(nsnap)
cd12 = np.zeros(nsnap)
ad12 = np.zeros(nsnap)
LS1 = np.zeros(nsnap)
LS2 = np.zeros(nsnap)
for i in range(nsnap):
    cs12[i] = get_angle(c_star1[i, :], c_star2[i, :])
    as12[i] = get_angle(a_star1[i, :], a_star2[i, :])
    cd12[i] = get_angle(c_dark1[i, :], c_dark2[i, :])
    ad12[i] = get_angle(a_dark1[i, :], a_dark2[i, :])
    LS1[i] = get_angle(L[i, :], spinVector1[i, :], acute=False)
    LS2[i] = get_angle(L[i, :], spinVector2[i, :], acute=False)

with open('{}/merger/orbit.dat'.format(args[0]), 'wb') as f:
    pickle.dump(rst, f)
with open('{}/merger/orbit.txt'.format(args[0]), 'w') as f:
    f.write('{:8s} {:>7s} {:>7s} {:>6s} {:>6s} {:>9s}  {:>5s}  {:>5s}  {:>5s}  {:>5s} {:>5s}'
            ' {:>5s} {:>5s} {:>5s} {:>5s} {:>5s}\n'
            .format('snapID', 'P1ID', 'P2ID', 'MRs', 'MRd', 'Lnorm', 'L1', 'L2', 'L3',
                    'V1', 'V2', 'V3', 'X1', 'X2', 'X3', 'r' ))
    for i in range(nsnap):
        f.write('{:8s} {:>7s} {:>7s} {:6.1f} {:6.1f} {:9.2e} [{:5.2f}, {:5.2f}, {:5.2f}]'
                ' {:5.0f} {:5.0f} {:5.0f} {:5.0f} {:5.0f} {:5.0f} {:5.0f}\n'
                .format(flist1[i, 0], flist1[i, 1][7:], flist2[i, 1][7:], mratio_star[i],
                        mratio_dark[i], Lnorm[i], L[i, 0], L[i, 1], L[i, 2], relative_vel[i, 0],
                        relative_vel[i, 2], relative_vel[i, 1], relative_pos[i, 0],
                        relative_pos[i, 1], relative_pos[i, 2], r[i]))

with open('{}/merger/plane.txt'.format(args[0]), 'w') as f:
    f.write('{:16s}  {:>5s}  {:>5s}  {:>5s}   {:>5s}  {:>5s}  {:>5s}   {:>5s}'
            '  {:>5s}  {:>5s}   {:>5s}  {:>5s}  {:>5s}  {:>5s} {:>5s} {:>5s}\n'
            .format('snapID', 'P1', 'P2', 'P3', 'L1', 'L2', 'L3',
                    'S11', 'S12', 'S13', 'S21', 'S22', 'S23', 'r', 'LS1', 'LS2'))
    for i in range(nsnap-1):
        temID = '{}-{}'.format(flist1[i, 0], flist1[i+1, 0])
        f.write('{:16s} [{:5.2f}, {:5.2f}, {:5.2f}] [{:5.2f}, {:5.2f}, {:5.2f}] '
                '[{:5.2f}, {:5.2f}, {:5.2f}] [{:5.2f}, {:5.2f}, {:5.2f}] {:5.0f} '
                '{:5.0f} {:5.0f}\n'
                .format(temID, orbit_plane[i, 0], orbit_plane[i, 1],
                        orbit_plane[i, 2], L[i, 0], L[i, 1], L[i, 2],
                        spinVector1[i, 0], spinVector1[i, 1], spinVector1[i, 2],
                        spinVector2[i, 0], spinVector2[i, 1], spinVector2[i, 2],
                        r[i], LS1[i], LS2[i]))

with open('{}/merger/angle_star.txt'.format(args[0]), 'w') as f:
    f.write('{:8s} {:>7s} {:>7s}  {:>5s}  {:>5s}  {:>5s}   {:>5s}  {:>5s}  {:>5s}  {:>5s}  '
            '{:>5s}  {:>5s}  {:>5s}   {:>5s}  {:>5s}  {:>5s}  {:>5s}\n'
            .format('snapID', 'P1ID', 'P2ID', 'cs1-1', 'cs1-2', 'cs1-3', 'cs2-1', 'cs2-2',
                    'cs2-3', 'Acs12', 'as1-1', 'as1-2', 'as1-3', 'as2-1', 'as2-2', 'as2-3',
                    'Aas12'))
    for i in range(nsnap):
        f.write('{:8s} {:>7s} {:>7s} [{:5.2f}, {:5.2f}, {:5.2f}] [{:5.2f}, {:5.2f}, {:5.2f}] '
                '{:5.0f} [{:5.2f}, {:5.2f}, {:5.2f}] [{:5.2f}, {:5.2f}, {:5.2f}] {:5.0f}\n'
                .format(flist1[i, 0], flist1[i, 1][7:], flist2[i, 1][7:], c_star1[i, 0],
                        c_star1[i, 1], c_star1[i, 2], c_star2[i, 0], c_star2[i, 1],
                        c_star2[i, 2], cs12[i], a_star1[i, 0], a_star1[i, 1], a_star1[i, 2],
                        a_star2[i, 0], a_star2[i, 1], a_star2[i, 2], as12[i]))

with open('{}/merger/angle_dark.txt'.format(args[0]), 'w') as f:
    f.write('{:8s} {:>7s} {:>7s}  {:>5s}  {:>5s}  {:>5s}   {:>5s}  {:>5s}  {:>5s}  {:>5s}  '
            '{:>5s}  {:>5s}  {:>5s}   {:>5s}  {:>5s}  {:>5s}  {:>5s}\n'
            .format('snapID', 'P1ID', 'P2ID', 'cd1-1', 'cd1-2', 'cd1-3', 'cd2-1', 'cd2-2',
                    'cd2-3', 'Acd12', 'ad1-1', 'ad1-2', 'ad1-3', 'ad2-1', 'ad2-2', 'ad2-3',
                    'Aad12'))
    for i in range(nsnap):
        f.write('{:8s} {:>7s} {:>7s} [{:5.2f}, {:5.2f}, {:5.2f}] [{:5.2f}, {:5.2f}, {:5.2f}] '
                '{:5.0f} [{:5.2f}, {:5.2f}, {:5.2f}] [{:5.2f}, {:5.2f}, {:5.2f}] {:5.0f}\n'
                .format(flist1[i, 0], flist1[i, 1][7:], flist2[i, 1][7:], c_dark1[i, 0],
                        c_dark1[i, 1], c_dark1[i, 2], c_dark2[i, 0], c_dark2[i, 1],
                        c_dark2[i, 2], cd12[i], a_dark1[i, 0], a_dark1[i, 1], a_dark1[i, 2],
                        a_dark2[i, 0], a_dark2[i, 1], a_dark2[i, 2], ad12[i]))
