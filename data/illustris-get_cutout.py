#!/usr/bin/env python
from optparse import OptionParser
import util_illustris as ui
import sys
import os
import pickle


parser = OptionParser()
parser.add_option('-s', action='store', type='int', dest='snapNum',
                  default=-1, help='snapshot number')
parser.add_option('-i', action='store', type='int', dest='subhaloID',
                  default=-1, help='subhalo ID')
(options, args) = parser.parse_args()

if options.snapNum == -1 or options.subhaloID == -1:
    print 'Error - please provide snapNum and subhaloID'
    sys.exit(1)

snapNum = options.snapNum
subID = options.subhaloID
os.system('mkdir -p snap{}-subhalo{}'.format(snapNum, subID))
url = ('http://www.illustris-project.org/api/Illustris-1/'
       'snapshots/{}/subhalos/{}/'.format(snapNum, subID))
cutout_url = ui.http_get(url)['cutouts']['subhalo']
info_url = ui.http_get(url)['meta']['info']
info = ui.http_get(info_url)
info_dict = ui.http_info2dict(info)
with open('snap{}-subhalo{}/info.dat'.format(snapNum, subID), 'wb') as f:
    pickle.dump(info_dict, f)

saved_filename = ui.http_get(cutout_url)
os.system('mv {} snap{}-subhalo{}/cutout.hdf5'
          .format(saved_filename, snapNum, subID))
print 'Download snap{}-subhalo{} success!'.format(snapNum, subID)
