#!/usr/bin/env python
import numpy as np
import os
from optparse import OptionParser
import sys


parser = OptionParser()
(options, args) = parser.parse_args()
if len(args) != 1:
    print 'Error - please provide a folder name'
    sys.exit(1)

flist = np.genfromtxt('{}/select.list'.format(args[0]), dtype='S30')
for i in range(flist.shape[0]):
    name = '{}-{}'.format(flist[i, 0], flist[i, 1])
    print name
    os.system('illustris-get_shape.py {}/{}'.format(args[0], name))
    os.system('illustris-plot_shape_single.py {}/{}'.format(args[0], name))
