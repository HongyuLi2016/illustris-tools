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
url = ('http://www.illustris-project.org/api/Illustris-1/'
       'snapshots/{}/subhalos/{}/'.format(snapNum, subID))
subhalo = ui.http_get(url)
url_halo = subhalo['related']['parent_halo']
url_halo_cutout = subhalo['cutouts']['parent_halo']
url_info = ui.http_get(url_halo)['meta']['info']
info = ui.http_get(url_info)
info_dict = ui.http_info2dict(info)

os.system('mkdir -p snap{}-halo{}'
          .format(info_dict['SnapshotNumber'], info_dict['InfoID']))
with open('snap{}-halo{}/info.dat'.format(info_dict['SnapshotNumber'],
                                          info_dict['InfoID']), 'wb') as f:
    pickle.dump(info_dict, f)
saved_filename = ui.http_get(url_halo_cutout)
os.system('mv {} snap{}-halo{}/cutout.hdf5'
          .format(saved_filename, info_dict['SnapshotNumber'],
                  info_dict['InfoID']))
print 'Download snap{}-halo{} success!'.format(info_dict['SnapshotNumber'],
                                               info_dict['InfoID'])
