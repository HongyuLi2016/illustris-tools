#!/usr/bin/env python
from optparse import OptionParser
import util_illustris as ui
import sys
import os
import numpy as np
import pickle


parser = OptionParser()
parser.add_option('-s', action='store', type='float', dest='size',
                  default=20.0, help='inclination')
parser.add_option('-r', action='store_true', dest='rotate',
                  default=False, help='rotate to principle axis coordinates')
parser.add_option('-l', action='store_true', dest='los',
                  default=False, help='make los velocity plot')
parser.add_option('-c', action='store_true', dest='cmark',
                  default=False, help='plot a marker for center')
(options, args) = parser.parse_args()
if len(args) != 1:
    print 'Error - please provide a folder name'
    sys.exit(1)

snapNum, subhaloID = ui.extract_id(args[0])
redshift = ui.snap2z(snapNum)
os.system('mkdir -p {}/info/'.format(args[0]))
size = options.size
box = np.array([[-size, size],
                [-size, size],
                [-size, size]]).T

cutout = ui.cutout('{}/cutout.hdf5'.format(args[0]))
'''
use info.dat
'''
with open('{}/info.dat'.format(args[0])) as f:
    info = pickle.load(f)
xcenter_star = np.array(info['Subhalo']['SubhaloPos'])
xcenter_dark = xcenter_star
hmr_star = info['Subhalo']['SubhaloHalfmassRadType'][4]
hmr_dark = info['Subhalo']['SubhaloHalfmassRadType'][1]
scale = 1.0

vpart_star = cutout.vpart_star
xpart_star = cutout.xpart_star
mpart_star = cutout.mass_star
# xcenter_star = ui.find_center(xpart_star, mpart=mpart_star, percent=10)

vpart_dark = cutout.vpart_dark
xpart_dark = cutout.xpart_dark
mpart_dark = cutout.mass_dark

# xcenter_dark = ui.find_center(xpart_dark, mpart=mpart_dark, percent=10)
center_offset = xcenter_star - xcenter_dark

vcenter_star = ui.vel_center(xpart_star-xcenter_star, vpart_star,
                             mpart=mpart_star, R=hmr_star)
vcenter_dark = ui.vel_center(xpart_dark-xcenter_dark, vpart_dark,
                             mpart=mpart_dark, R=hmr_dark)
vsys_offset = vcenter_star - vcenter_dark

# subtract center position and system velocity (position are based on xstar)
xpart_star -= xcenter_star
xpart_dark -= xcenter_star
vpart_star -= vcenter_star
vpart_dark -= vcenter_dark
# mark a potision on the image
xmark_star = np.array([[18612.98,  68966.09,  3290.10],
                       [18624.92,  68963.44,  3273.73]])
xmark_star = np.array([52143.1,  36207.2,  24010.1])
xmark_star -= xcenter_star
if not options.cmark:
    xmark_star = None

i_projection = [(0, 1, 0), (0, 2, 1), (1, 2, 2)]
if options.rotate:
    ba, ca, angle, Tiv = \
        ui.get_shape(xpart_star, mpart_star, Rb=scale*hmr_star, decrease=True)
    xpart_star = np.dot(Tiv, xpart_star.T).T
    vpart_star = np.dot(Tiv, vpart_star.T).T
    xpart_dark = np.dot(Tiv, xpart_dark.T).T
    vpart_dark = np.dot(Tiv, vpart_dark.T).T


fig_star, axes_star = ui.cutout_vel_vector(xpart_star, vpart_star, mpart_star,
                                           alpha=0.7, headlength=6.0, box=box,
                                           ifu_bins=50, xmark=xmark_star, magrange=7)
fig_star.savefig('{}/info/vel_vector_star.png'.format(args[0]), dpi=500)

fig_dark, axes_dark = ui.cutout_vel_vector(xpart_dark, vpart_dark, mpart_dark,
                                           alpha=0.7, headlength=6.0, box=box,
                                           ifu_bins=50, magrange=3,
                                           xmark=xmark_star)
fig_dark.savefig('{}/info/vel_vector_dark.png'.format(args[0]), dpi=500)


if options.los:
    fig_star, axes_star = ui.cutout_vel_los(xpart_star, vpart_star, mpart_star,
                                            box=box, linewidths=0.5,
                                            colors='c', xmark=xmark_star)
    fig_star.savefig('{}/info/vel_los_star.png'.format(args[0]), dpi=500)

    fig_dark, axes_dark = ui.cutout_vel_los(xpart_dark, vpart_dark, mpart_dark,
                                            box=box, linewidths=0.3,
                                            colors='c', xmark=xmark_star)
    fig_dark.savefig('{}/info/vel_los_dark.png'.format(args[0]), dpi=500)

with open('{}/info/cutoutInfo.txt'.format(args[0]), 'w') as ff:
    ff.write('SnapNum: {}\n'.format(snapNum))
    ff.write('Redshift: {:.4f}\n'.format(redshift))
    ff.write('SubhaloID: {}\n'.format(subhaloID))
    ff.write('Dark particle number: {:.4e}\n'.format(cutout.numDark))
    ff.write('Star particle number: {:.4e}\n'.format(cutout.numStar))
    ff.write('Dark mass (all particals): {:.3f}\n'.format(cutout.Mdark))
    ff.write('Star mass (all particals): {:.3f}\n'.format(cutout.Mstar))
    ff.write('Dark matter center: {:.2f}  {:.2f}  {:.2f}\n'
             .format(xcenter_dark[0], xcenter_dark[1], xcenter_dark[2]))
    ff.write('Star center: {:.2f}  {:.2f}  {:.2f}\n'
             .format(xcenter_star[0], xcenter_star[1], xcenter_star[2]))
    ff.write('Center offset: {:.2f}  {:.2f}  {:.2f}  {:.2f}\n'
             .format(center_offset[0], center_offset[1], center_offset[2],
                     np.sum(center_offset**2)**0.5))
    ff.write('Dark matter vsys: {:.2f}  {:.2f}  {:.2f}\n'
             .format(vcenter_dark[0], vcenter_dark[1], vcenter_dark[2]))
    ff.write('Star center vsys: {:.2f}  {:.2f}  {:.2f}\n'
             .format(vcenter_star[0], vcenter_star[1], vcenter_star[2]))
    ff.write('Vsys offset: {:.2f}  {:.2f}  {:.2f}  {:.2f}\n'
             .format(vsys_offset[0], vsys_offset[1], vsys_offset[2],
                     np.sum(vsys_offset**2)**0.5))
