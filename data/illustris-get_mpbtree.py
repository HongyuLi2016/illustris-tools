#!/usr/bin/env python
from optparse import OptionParser
import util_illustris as ui
import sys
import os


parser = OptionParser()
(options, args) = parser.parse_args()
if len(args) != 1:
    print 'Error - please provide a folder name'
    sys.exit(1)
os.system('mkdir -p {}'.format(args[0]))
ID = args[0][7:]
url = ('http://www.illustris-project.org/api/Illustris-1/'
       'snapshots/135/subhalos/{}/sublink/full.hdf5'.format(ID))
saved_filename = ui.http_get(url)
os.system('mv {} {}/sublink_tree.hdf5'.format(saved_filename, args[0]))
print 'Download {} success!'.format(args[0])
