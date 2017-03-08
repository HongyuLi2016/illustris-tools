#!/usr/bin/env python
'''
select progenitors to be checked for a given merger node
'''
import numpy as np
from illustris_tree_basic import tree_basic
from optparse import OptionParser
import os


if __name__ == '__main__':
    parser = OptionParser()
    (options, args) = parser.parse_args()
    if len(args) != 1:
        print 'Error - please provide a folder name'
        sys.exit(1)

    os.system('mkdir -p {}/merger/progenitor1/basic_plots'.format(args[0]))
    os.system('mkdir -p {}/merger/progenitor2/basic_plots'.format(args[0]))

    node = np.genfromtxt('{}/node-merger'.format(args[0]), dtype=int)
    snap1 = node[0, 0]
    ID1 = node[0, 1]
    snap2 = node[1, 0]
    ID2 = node[1, 1]

    lhy1 = tree_basic(args[0], rootSnap=snap1,
                      rootSubfindID=ID1)
    lhy1.dump(outpath='{}/merger/progenitor1'.format(args[0]))
    lhy1.dumpTxt(outpath='{}/merger/progenitor1'.format(args[0]))
    lhy1.plot(outpath='{}/merger/progenitor1/basic_plots'.format(args[0]))

    lhy2 = tree_basic(args[0], rootSnap=snap2,
                      rootSubfindID=ID2)
    lhy2.dump(outpath='{}/merger/progenitor2'.format(args[0]))
    lhy2.dumpTxt(outpath='{}/merger/progenitor2'.format(args[0]))
    lhy2.plot(outpath='{}/merger/progenitor2/basic_plots'.format(args[0]))

    flist1 = np.genfromtxt('{}/merger/progenitor1/mpbTree.txt'.format(args[0]),
                           usecols=[0, 2], dtype='S30')
    flist2 = np.genfromtxt('{}/merger/progenitor2/mpbTree.txt'.format(args[0]),
                           usecols=[0, 2], dtype='S30')
    in1 = np.in1d(flist1[:, 0], flist2[:, 0])
    in2 = np.in1d(flist2[:, 0], flist1[:, 0])
    common_list1 = flist1[in1, :]
    common_list2 = flist2[in2, :]
    with open('{}/merger/progenitor1/select.list'.format(args[0]), 'w') as f:
        for i in range(5):
            f.write('{} {}\n'.format(common_list1[i, 0], common_list1[i, 1]))
    with open('{}/merger/progenitor2/select.list'.format(args[0]), 'w') as f:
        for i in range(5):
            f.write('{} {}\n'.format(common_list2[i, 0], common_list2[i, 1]))
