#!/usr/bin/env python
from optparse import OptionParser
import util_illustris as ui
import sys
import os


parser = OptionParser()
parser.add_option('-s', action='store', type='int', dest='snapNum',
                  default=-1, help='snapshot number')
parser.add_option('-i', action='store', type='int', dest='subhaloID',
                  default=-1, help='subhalo ID')
(options, args) = parser.parse_args()
if len(args) != 1:
    print 'Error - please provide a folder name'
    sys.exit(1)

snapNum = options.snapNum
subID = options.subhaloID
os.system('mkdir -p {}/snap{}-subhalo{}'.format(args[0], snapNum, subID))
url = ('http://www.illustris-project.org/api/Illustris-1/'
       'snapshots/{}/subhalos/{}/cutout.hdf5'.format(snapNum, subID))
saved_filename = ui.http_get(url)
os.system('mv {} {}/snap{}-subhalo{}/cutout.hdf5'
          .format(saved_filename, args[0], snapNum, subID))
print 'Download {}/snap{}-subhalo{} success!'.format(args[0], snapNum, subID)
