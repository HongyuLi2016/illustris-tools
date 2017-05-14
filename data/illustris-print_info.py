#!/usr/bin/env python
from optparse import OptionParser
import sys
import pickle

parser = OptionParser()
(options, args) = parser.parse_args()
if len(args) != 1:
    print 'Error - please provide snapNum and subhaloID'
    sys.exit(1)
with open('{}/info.dat'.format(args[0]), 'rb') as f:
    info = pickle.load(f)
for key in info['Subhalo'].keys():
    print key, info['Subhalo'][key]

with open('{}/haloInfo.dat'.format(args[0]), 'rb') as f:
    info = pickle.load(f)
for key in info['Group'].keys():
    print key, info['Group'][key]

