#!/usr/bin/env python
from illustris_tree_basic import tree_basic
from optparse import OptionParser
import sys
import os


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-s', action='store', type='int', dest='snapNum',
                      default=None, help='snapshot number')
    parser.add_option('-i', action='store', type='int', dest='subhaloID',
                      default=None, help='subhalo ID')

    (options, args) = parser.parse_args()
    if len(args) != 1:
        print 'Error - please provide a folder name'
        sys.exit(1)
    os.system('mkdir -p {}/basic_plots'.format(args[0]))
    lhy = tree_basic(args[0], rootSnap=options.snapNum,
                     rootSubfindID=options.subhaloID)
    lhy.dump(outpath=args[0])
    lhy.dumpTxt(outpath=args[0])
    lhy.plot(outpath='{}/basic_plots'.format(args[0]))
