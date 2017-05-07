#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pickle
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-l', action='store', type='float', dest='lim',
                  default=np.nan, help='box size')
parser.add_option('-n', action='store', type='int', dest='n',
                  default=-1, help='plot second progenitor')
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
S1 = orbit['spinVector1']
S2 = orbit['spinVector2']
try:
    good = np.load('{}/merger/good.npy'.format(args[0]))
except:
    print 'No good file provided'
    good = np.ones(len(Lnorm)-1, dtype=bool)
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
length = 0.2 * lim
ax.quiver([0], [0], [0], S1[1:, 0][good], S1[1:, 1][good],
          S1[1:, 2][good], pivot='tail', length=length, color='y')
n = options.n
ax.quiver([r[n, 0]], [r[n, 1]], [r[n, 2]], S2[1:, 0][good], S2[1:, 1][good],
          S2[1:, 2][good], pivot='tail', length=length, color='c')
ax.quiver(r[:, 0], r[:, 1], r[:, 2], v[:, 0], v[:, 1],
          v[:, 2], pivot='tail', length=length)
ax.set_aspect(1)
plt.show()
