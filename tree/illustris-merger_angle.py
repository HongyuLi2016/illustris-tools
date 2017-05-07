#!/usr/bin/env python
import numpy as np
import itertools
from optparse import OptionParser
import os
import pickle


def get_angle(vect1, vect2):
    '''
    Calculate the angle between 2 vectors
    '''
    vect1 = vect1 / np.linalg.norm(vect1)
    vect2 = vect2 / np.linalg.norm(vect2)
    cosi = np.dot(vect1, vect2)
    angle = np.degrees(np.arccos(cosi))
    # if angle >= 90.0:
    #     angle = 180 - angle
    return angle


parser = OptionParser()
(options, args) = parser.parse_args()
if len(args) != 1:
    print 'Error - please provide a folder name'
    exit(1)
with open('{}/merger/orbit.dat'.format(args[0]), 'rb') as f:
    rst = pickle.load(f)
as1 = rst['a_star1']
as2 = rst['a_star2']
ad1 = rst['a_dark1']
ad2 = rst['a_dark2']
cs1 = rst['c_star1']
cs2 = rst['c_star2']
cd1 = rst['c_dark1']
cd2 = rst['c_dark2']
S1 = rst['spinVector1']
S2 = rst['spinVector2']
L = rst['L']
vectors = {'as1': as1, 'as2': as2, 'ad1': ad1, 'ad2': ad2, 'cs1': cs1,
           'cs2': cs2, 'cd1': cd1, 'cd2': cd2, 'L': L, 'S1': S1, 'S2':S2}
names = itertools.combinations(vectors.keys(), 2)
angles = {}
for name in names:
    comb_names = '{}-{}'.format(name[0], name[1])
    angles[comb_names] = \
        [get_angle(vectors[name[0]][i, :], vectors[name[1]][i, :])
         for i in range(L.shape[0])]
with open('{}/merger/angles.dat'.format(args[0]), 'wb') as f:
    pickle.dump(angles, f)
