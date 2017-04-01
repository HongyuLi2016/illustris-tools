#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pickle
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-a', action='store_false', dest='c',
                  default=True, help='plot a axis')
parser.add_option('-d', action='store_false', dest='star',
                   default=True, help='plot dark matter')
parser.add_option('-l', action='store', type='float', dest='lim',
                   default=np.nan, help='box size')
(options, args) = parser.parse_args()
if len(args) != 1:
    print 'Error - please provide a folder name'
    exit(1)


with open('{}/merger/orbit.dat'.format(args[0])) as f:
    orbit = pickle.load(f)
r = orbit['relative_pos']
v = orbit['relative_vel']
Lnorm = orbit['Lnorm']
a_star1 = orbit['a_star1']
a_star2 = orbit['a_star2']
c_star1 = orbit['c_star1']
c_star2 = orbit['c_star2']
a_dark1 = orbit['a_dark1']
a_dark2 = orbit['a_dark2']
c_dark1 = orbit['c_dark1']
c_dark2 = orbit['c_dark2']
fig = plt.figure(figsize=(15, 15))
ax = fig.gca(projection='3d')
ax.plot([0], [0], [0], 'or')
ax.plot(r[:, 0], r[:, 1], r[:, 2], 'o-k')
ax.plot([r[0, 0]], [r[0, 1]], [r[0, 2]], 'o-g')
ax.plot([r[-1, 0]], [r[-1, 1]], [r[-1, 2]], 'o-y')
if not np.isnan(options.lim):
    lim = options.lim
else:
    lim = np.max(np.abs(r)) * 0.9
ax.set_xlim([-lim, lim])
ax.set_ylim([-lim, lim])
ax.set_zlim([-lim, lim])
length = 0.3 * lim
print 'blue: progenitor 1  green: progenitor 2'
if options.star:
    if options.c:
        ax.quiver(r[:, 0], r[:, 1], r[:, 2], c_star1[:, 0], c_star1[:, 1],
                  c_star1[:, 2], pivot='tail', length=length, color='b')
        ax.quiver(r[:, 0], r[:, 1], r[:, 2], c_star2[:, 0], c_star2[:, 1],
                  c_star2[:, 2], pivot='tail', length=length, color='g')
    else:
        ax.quiver(r[:, 0], r[:, 1], r[:, 2], a_star1[:, 0], a_star1[:, 1],
                  a_star1[:, 2], color='b', pivot='tail', length=length)
        ax.quiver(r[:, 0], r[:, 1], r[:, 2], a_star2[:, 0], a_star2[:, 1],
                  a_star2[:, 2], color='g', pivot='tail', length=length)
else:
    print 'Plot dark matter axis'
    if options.c:
        ax.quiver(r[:, 0], r[:, 1], r[:, 2], c_dark1[:, 0], c_dark1[:, 1],
                  c_dark1[:, 2], pivot='tail', length=length, color='b')
        ax.quiver(r[:, 0], r[:, 1], r[:, 2], c_dark2[:, 0], c_dark2[:, 1],
                  c_dark2[:, 2], pivot='tail', length=length, color='g')
    else:
        ax.quiver(r[:, 0], r[:, 1], r[:, 2], a_dark1[:, 0], a_dark1[:, 1],
                  a_dark1[:, 2], color='b', pivot='tail', length=length)
        ax.quiver(r[:, 0], r[:, 1], r[:, 2], a_dark2[:, 0], a_dark2[:, 1],
                  a_dark2[:, 2], color='g', pivot='tail', length=length)
ax.set_aspect(1)
plt.show()
