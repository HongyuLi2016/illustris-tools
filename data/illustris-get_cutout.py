#!/usr/bin/env python
'''
Download a particle cut out file for a given ID and snapshot
Small: instead of downloading all the field, just download the
stellar and dark matter particle position, velocity, mass and
metallicity
'''
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

parentHaloUrl = ui.http_get(url)['related']['parent_halo']
parentHaloInfoUrl = ui.http_get(parentHaloUrl)['meta']['info']
haloInfo = ui.http_get(parentHaloInfoUrl)
haloInfo_dict = ui.http_info2dict(haloInfo)

with open('snap{}-subhalo{}/info.dat'.format(snapNum, subID), 'wb') as f:
    pickle.dump(info_dict, f)
with open('snap{}-subhalo{}/haloInfo.dat'.format(snapNum, subID), 'wb') as f:
    pickle.dump(haloInfo_dict, f)

params = {'stars': 'Coordinates,Velocities,Masses,GFM_Metallicity,'
          'GFM_StellarPhotometrics',
          'dm': 'Coordinates,Velocities'}
success = False
try_num = 0
f = open('snap{}-subhalo{}/download.log'.format(snapNum, subID), 'w')
while not success:
    try:
        if try_num > 5:
            break
        try_num += 1
        saved_filename = ui.http_get(cutout_url, params)
        success = True
        os.system('mv {} snap{}-subhalo{}/cutout.hdf5'
                  .format(saved_filename, snapNum, subID))
        print>>f, 'Download snap{}-subhalo{} success!'.format(snapNum, subID)
    except:
        print>>f, 'Try {} failed!'.format(try_num)
f.close()
