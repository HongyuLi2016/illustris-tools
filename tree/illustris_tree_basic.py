#!/usr/bin/env python
# import numpy as np
import h5py
import util_illustris as ui
import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import colors, colorbar
import warnings
from JAM.utils import util_fig
util_fig.ticks_font.set_size(15)
warnings.simplefilter("ignore")


def normRadius(R, R11=0.3, R9=0.05):
    Rnorm = (R11-R9)/(11.0-9.0) * (R - 9.0) + R9
    return Rnorm


def nodeCopy(nodeDest, nodeOrigin):
    if set(nodeDest.keys()) != set(nodeOrigin.keys()):
        raise Warning('Keys not equal in a dictionary copy')
    for key in nodeDest:
        nodeDest[key] = nodeOrigin[key]


def nodeNew(node):
    new = {}
    for key in node:
        new[key] = node[key]
    return new


class tree_basic_original_tree:
    '''
    trace the original mbp tree given by sublink
    '''
    def __init__(self, path, rootSnap=None, rootSubfindID=None):
        self.path = path
        self.f = h5py.File('{}/sublink_tree.hdf5'.format(path), 'r')

        self.FirstProgenitorID = self.get('FirstProgenitorID')
        self.NextProgenitorID = self.get('NextProgenitorID')
        self.SubhaloID = self.get('SubhaloID')
        self.SnapNum = self.get('SnapNum')
        self.SubfindID = self.get('SubfindID')
        self.SubhaloCM = self.get('SubhaloCM')
        self.SubhaloHalfmassRad = self.get('SubhaloHalfmassRad')
        self.SubhaloHalfmassRadType = self.get('SubhaloHalfmassRadType')
        self.SubhaloLenType = self.get('SubhaloLenType')
        self.SubhaloMass = self.get('SubhaloMass')
        self.SubhaloMassType = self.get('SubhaloMassType')
        self.SubhaloMassInHalfRad = self.get('SubhaloMassInHalfRad')
        self.SubhaloMassInHalfRadType = self.get('SubhaloMassInHalfRadType')
        self.SubhaloPos = self.get('SubhaloPos')
        self.SubhaloSFR = self.get('SubhaloSFR')
        self.SubhaloSFRinHalfRad = self.get('SubhaloSFRinHalfRad')
        self.SubhaloVel = self.get('SubhaloVel')
        self.Mass = self.get('Mass')
        self.MassHistory = self.get('MassHistory')
        # create the first node
        if rootSnap is None and rootSubfindID is not None:
            raise ValueError('both rootSnap and'
                             ' rootSubfindID must be provided!')
        if rootSnap is not None and rootSubfindID is None:
            raise ValueError('both rootSnap and'
                             ' rootSubfindID must be provided!')
        if rootSnap is None:
            rootSubfindID = self.SubfindID[0]
            rootSnap = self.SnapNum[0]
        iRoot = (self.SnapNum == rootSnap) * (self.SubfindID == rootSubfindID)
        if iRoot.sum() < 1:
            raise RuntimeError('Subhalo not found')
        if iRoot.sum() > 1:
            raise RuntimeError('More than 1 subhalos are found')
        rootID = self.SubhaloID[iRoot][0]  # subhaloID for root node
        self.tree = self.create_node(rootID)  # create root node
        node = self.tree
        nodeID = rootID
        while True:
            nodeIndex = nodeID - self.SubhaloID[0]
            if self.FirstProgenitorID[nodeIndex] == -1:
                break
            if self.NextProgenitorID[nodeIndex] == -1:
                node['NextProgenitor'] = None  # create next progenitor
                firstNodeID = self.FirstProgenitorID[nodeIndex]
            else:
                # create next progenitor
                nextNodeID = self.NextProgenitorID[nodeIndex]
                node['NextProgenitor'] = self.create_node(nextNodeID)
                # Select firstNodeID
                # choose between the 'first progenitor' and the 'next progenitor
                # of the fisrt progenitor'
                # (subhalo with larger Mstar as the first progenitor)
                firstNodeID1 = self.FirstProgenitorID[nodeIndex]
                firstIndex1 = firstNodeID1 - self.SubhaloID[0]
                firstNodeID2 = self.NextProgenitorID[firstIndex1]
                firstIndex2 = firstNodeID2 - self.SubhaloID[0]
                if firstIndex2 < 0:
                    firstNodeID = firstNodeID1
                else:
                    # firstNode1Mstar = self.SubhaloMassType[firstIndex1][4]
                    # firstNode2Mstar = self.SubhaloMassType[firstIndex2][4]
                    firstNode1Mstar = self.MassHistory[firstIndex1]
                    firstNode2Mstar = self.MassHistory[firstIndex2]
                    if firstNode1Mstar > firstNode2Mstar:
                        firstNodeID = firstNodeID1
                    else:
                        '''
                        print ('Waring - StellarMass First progenior < '
                               'StellarMass Next Progenitor for snapNum {}'
                               'Fsubhalo{} Nsubhalo{}'
                               .format(node['SnapNum'], node['SubfindID'],
                               Nextnode['SubfindID']))
                        '''
                        firstNodeID = firstNodeID2
            # create first progenitor node and link
            node['FirstProgenitor'] = self.create_node(firstNodeID)
            nodeID = firstNodeID
            node = node['FirstProgenitor']
        node['FirstProgenitor'] = None
        node['NextProgenitor'] = None

    def get(self, key):
        return self.f[key][:]

    def create_node(self, subhaloID):
        node = {}
        index = subhaloID - self.SubhaloID[0]
        if index < 0:
            raise RuntimeError('subhaloID too small')
        node['SnapNum'] = self.SnapNum[index]
        node['SubhaloID'] = subhaloID
        node['redshift'] = ui.snap2z(node['SnapNum'])
        node['SubfindID'] = self.SubfindID[index]
        # info from SUBFIND catalogue
        node['SubhaloCM'] = self.SubhaloCM[index]
        node['SubhaloHalfmassRad'] = self.SubhaloHalfmassRad[index]
        node['SubhaloHalfmassRadType'] = self.SubhaloHalfmassRadType[index]
        node['SubhaloLenType'] = self.SubhaloLenType[index]
        node['SubhaloMass'] = self.SubhaloMass[index]
        node['SubhaloMassType'] = self.SubhaloMassType[index]
        node['SubhaloMassInHalfRad'] = self.SubhaloMassInHalfRad[index]
        node['SubhaloMassInHalfRadType'] = self.SubhaloMassInHalfRadType[index]
        node['SubhaloPos'] = self.SubhaloPos[index]
        node['SubhaloSFR'] = self.SubhaloSFR[index]
        node['SubhaloSFRinHalfRad'] = self.SubhaloSFRinHalfRad[index]
        node['SubhaloVel'] = self.SubhaloVel[index]
        node['Mass'] = self.Mass[index]
        node['MassHistory'] = self.MassHistory[index]
        return node

    def dump(self, fname='mpbTree.dat', outpath='.'):
        with open('{}/{}'.format(outpath, fname), 'wb') as f:
            pickle.dump(self.tree, f)

    def dumpTxt(self, fname='mpbTree.txt', outpath='.'):
        node = self.tree
        with open('{}/{}'.format(outpath, fname), 'wb') as f:
            while node['FirstProgenitor'] is not None:
                mhalo = np.log10(node['SubhaloMassType'][1]*ui.massUnit)
                mstar = np.log10(node['SubhaloMassType'][4]*ui.massUnit)
                MassHistory = np.log10(node['MassHistory']*ui.massUnit)
                pos = node['SubhaloPos']
                # mass = np.log10(node['Mass']*ui.massUnit)
                if node['NextProgenitor'] is not None:
                    nextID = node['NextProgenitor']['SubfindID']
                    nextSnap = node['NextProgenitor']['SnapNum']
                    nextMhalo =\
                        np.log10(node['NextProgenitor']['SubhaloMassType'][1] *
                                 ui.massUnit)
                    nextMstar =\
                        np.log10(node['NextProgenitor']['SubhaloMassType'][4] *
                                 ui.massUnit)
                    nextMassHistory =\
                        np.log10(node['NextProgenitor']['MassHistory'] *
                                 ui.massUnit)
                    nextPos = node['NextProgenitor']['SubhaloPos']
                    nextstr = ('snap{:<3d}   subhalo{:<8d} {:7.2f} {:7.2f}'
                               '{:7.2f} [{:8.1f}, {:8.1f}, {:8.1f}]\n'
                               .format(nextSnap, nextID, nextMhalo, nextMstar,
                                       nextMassHistory, nextPos[0], nextPos[1],
                                       nextPos[2]))
                else:
                    nextstr = '\n'
                f.write('snap{:<3d} {:6.3f}   subhalo{:<8d} {:7.2f} {:7.2f}'
                        '{:7.2f} [{:8.1f}, {:8.1f}, {:8.1f}]   {}'
                        .format(node['SnapNum'], node['redshift'],
                                node['SubfindID'], mhalo, mstar, MassHistory,
                                pos[0], pos[1], pos[2], nextstr))
                node = node['FirstProgenitor']


class tree_basic:
    '''
    trace the stellar mass mbp tree used in prolate project
    '''
    def __init__(self, path, rootSnap=None, rootSubfindID=None):
        self.path = path
        self.f = h5py.File('{}/sublink_tree.hdf5'.format(path), 'r')

        self.DescendantID = self.get('DescendantID')
        self.FirstProgenitorID = self.get('FirstProgenitorID')
        self.NextProgenitorID = self.get('NextProgenitorID')
        self.SubhaloID = self.get('SubhaloID')
        self.SnapNum = self.get('SnapNum')
        self.SubfindID = self.get('SubfindID')
        self.SubhaloCM = self.get('SubhaloCM')
        self.SubhaloHalfmassRad = self.get('SubhaloHalfmassRad')
        self.SubhaloHalfmassRadType = self.get('SubhaloHalfmassRadType')
        self.SubhaloLenType = self.get('SubhaloLenType')
        self.SubhaloMass = self.get('SubhaloMass')
        self.SubhaloMassType = self.get('SubhaloMassType')
        self.Mstar = self.SubhaloMassType[:, 4]
        self.SubhaloMassInHalfRad = self.get('SubhaloMassInHalfRad')
        self.SubhaloMassInHalfRadType = self.get('SubhaloMassInHalfRadType')
        self.SubhaloMassInRad = self.get('SubhaloMassInRad')
        self.SubhaloMassInRadType = self.get('SubhaloMassInRadType')
        self.SubhaloPos = self.get('SubhaloPos')
        self.SubhaloSFR = self.get('SubhaloSFR')
        self.SubhaloSFRinHalfRad = self.get('SubhaloSFRinHalfRad')
        self.SubhaloVel = self.get('SubhaloVel')
        self.Mass = self.get('Mass')
        self.MassHistory = self.get('MassHistory')
        # create the first node
        if rootSnap is None and rootSubfindID is not None:
            raise ValueError('both rootSnap and'
                             ' rootSubfindID must be provided!')
        if rootSnap is not None and rootSubfindID is None:
            raise ValueError('both rootSnap and'
                             ' rootSubfindID must be provided!')
        if rootSnap is None:
            rootSubfindID = self.SubfindID[0]
            rootSnap = self.SnapNum[0]
        iRoot = (self.SnapNum == rootSnap) * (self.SubfindID == rootSubfindID)
        if iRoot.sum() < 1:
            raise RuntimeError('Subhalo not found')
        if iRoot.sum() > 1:
            raise RuntimeError('More than 1 subhalos are found')
        rootID = self.SubhaloID[iRoot][0]  # subhaloID for root node
        self.tree = self.create_node(rootID)  # create root node
        node = self.tree
        node['next'] = None
        nodeID = rootID
        while True:
            iProgenitor = self.DescendantID == nodeID
            if iProgenitor.sum() == 0:
                break
            Mstar = self.Mstar[iProgenitor]
            progenitorID = self.SubhaloID[iProgenitor]
            sortIndex = np.argsort(Mstar)[::-1]
            Mstar = Mstar[sortIndex]
            progenitorID = progenitorID[sortIndex]
            firstNodeID = progenitorID[0]
            Progenitors = self.create_node(firstNodeID)
            temp = Progenitors
            for i in range(len(sortIndex)-1):
                if Mstar[i+1]/Mstar[0] > 0.01:
                    nextNodeID = progenitorID[i+1]
                    temp['next'] = self.create_node(nextNodeID)
                    temp = temp['next']
                else:
                    break
            temp['next'] = None
            node['Progenitors'] = Progenitors
            nodeID = firstNodeID
            node = Progenitors
        node['Progenitors'] = None

    def get(self, key):
        return self.f[key][:]

    def create_node(self, subhaloID):
        node = {}
        index = subhaloID - self.SubhaloID[0]
        if index < 0:
            raise RuntimeError('subhaloID too small')
        node['SnapNum'] = self.SnapNum[index]
        node['SubhaloID'] = subhaloID
        node['redshift'] = ui.snap2z(node['SnapNum'])
        node['SubfindID'] = self.SubfindID[index]
        # info from SUBFIND catalogue
        node['SubhaloCM'] = self.SubhaloCM[index]
        node['SubhaloHalfmassRad'] = self.SubhaloHalfmassRad[index]
        node['SubhaloHalfmassRadType'] = self.SubhaloHalfmassRadType[index]
        node['SubhaloLenType'] = self.SubhaloLenType[index]
        node['SubhaloMass'] = self.SubhaloMass[index]
        node['SubhaloMassType'] = self.SubhaloMassType[index]
        node['SubhaloMassInHalfRad'] = self.SubhaloMassInHalfRad[index]
        node['SubhaloMassInHalfRadType'] = self.SubhaloMassInHalfRadType[index]
        node['SubhaloMassInRad'] = self.SubhaloMassInRad[index]
        node['SubhaloMassInRadType'] = self.SubhaloMassInRadType[index]
        node['SubhaloPos'] = self.SubhaloPos[index]
        node['SubhaloSFR'] = self.SubhaloSFR[index]
        node['SubhaloSFRinHalfRad'] = self.SubhaloSFRinHalfRad[index]
        node['SubhaloVel'] = self.SubhaloVel[index]
        node['Mass'] = self.Mass[index]
        node['MassHistory'] = self.MassHistory[index]
        return node

    def dump(self, fname='mpbTree.dat', outpath='.'):
        with open('{}/{}'.format(outpath, fname), 'wb') as f:
            pickle.dump(self.tree, f)

    def dumpTxt(self, fname='mpbTree.txt', outpath='.'):
        node = self.tree
        with open('{}/{}'.format(outpath, fname), 'wb') as f:
            while True:
                mhalo = np.log10(node['SubhaloMassType'][1]*ui.massUnit)
                mstar = np.log10(node['SubhaloMassType'][4]*ui.massUnit)
                MassHistory = np.log10(node['MassHistory']*ui.massUnit)
                pos = node['SubhaloPos']
                # mass = np.log10(node['Mass']*ui.massUnit)
                if node['next'] is not None:
                    nextID = node['next']['SubfindID']
                    nextSnap = node['next']['SnapNum']
                    nextMhalo =\
                        np.log10(node['next']['SubhaloMassType'][1] *
                                 ui.massUnit)
                    nextMstar =\
                        np.log10(node['next']['SubhaloMassType'][4] *
                                 ui.massUnit)
                    nextMassHistory =\
                        np.log10(node['next']['MassHistory'] *
                                 ui.massUnit)
                    nextPos = node['next']['SubhaloPos']
                    nextstr = ('snap{:<3d}   subhalo{:<8d} {:7.2f} {:7.2f}'
                               '{:7.2f} [{:8.1f}, {:8.1f}, {:8.1f}]\n'
                               .format(nextSnap, nextID, nextMhalo, nextMstar,
                                       nextMassHistory, nextPos[0], nextPos[1],
                                       nextPos[2]))
                else:
                    nextstr = '\n'
                f.write('snap{:<3d} {:6.3f}   subhalo{:<8d} {:7.2f} {:7.2f}'
                        '{:7.2f} [{:8.1f}, {:8.1f}, {:8.1f}]   {}'
                        .format(node['SnapNum'], node['redshift'],
                                node['SubfindID'], mhalo, mstar, MassHistory,
                                pos[0], pos[1], pos[2], nextstr))
                if node['Progenitors'] is None:
                    break
                else:
                    node = node['Progenitors']

    def plot(self, outpath='.', start=60):
        node = self.tree
        snapNum = []
        mhalo = []
        mstar = []
        fdm = []
        fgas = []
        snapNumNext = []
        # mhaloNext = []
        mstarNext = []
        # fdmNext = []
        fgasNext = []
        while True:
            snapNum.append(node['SnapNum'])
            mhalo.append(np.log10(node['SubhaloMassType'][1]*ui.massUnit))
            mstar.append(np.log10(node['SubhaloMassType'][4]*ui.massUnit))
            mgasRad = node['SubhaloMassInRadType'][0]
            mdarkRad = node['SubhaloMassInRadType'][1]
            mstarRad = node['SubhaloMassInRadType'][4]
            fdm.append(mdarkRad / (mgasRad + mdarkRad + mstarRad))
            fgas.append(mgasRad / (mgasRad + mdarkRad + mstarRad))
            if node['next'] is not None:
                snapNumNext.append(node['next']['SnapNum'])
                mstarNext.append(np.log10(node['next']['SubhaloMassType'][4] *
                                 ui.massUnit))
                mgasRad = node['next']['SubhaloMassInRadType'][0]
                mdarkRad = node['next']['SubhaloMassInRadType'][1]
                mstarRad = node['next']['SubhaloMassInRadType'][4]
                fgasNext.append(mgasRad / (mgasRad + mdarkRad + mstarRad))
            else:
                snapNumNext.append(np.nan)
                mstarNext.append(np.nan)
                fgasNext.append(np.nan)
            if node['Progenitors'] is None:
                break
            else:
                node = node['Progenitors']
        snapNum = np.array(snapNum)
        mhalo = np.array(mhalo)
        mstar = np.array(mstar)
        fdm = np.array(fdm)
        fgas = np.array(fgas)
        snapNumNext = np.array(snapNumNext)
        mstarNext = np.array(mstarNext)
        fgasNext = np.array(fgasNext)
        plt.clf()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(snapNum, mhalo, 'k', label='dark')
        ax.plot(snapNum, mstar, 'b', label='star')
        ax.legend(loc='lower right')
        ax.set_xlabel('snapshot')
        ax.set_ylabel(r'Mass [$\rm \log \ M_\odot$]')
        fig.savefig('{}/MG.png'.format(outpath))
        plt.clf()
        ax = fig.add_subplot(111)
        ax.plot(snapNum, fdm, 'k', label='fdm(<2Hmr)')
        ax.plot(snapNum, fgas, 'g', label='gas(<2Hmr)')
        ax.legend(loc='lower right')
        ax.set_xlabel('snapshot')
        ax.set_ylabel(r'Mass [$\rm \log \ M_\odot$]')
        fig.savefig('{}/gas-fdm.png'.format(outpath))
        ii = snapNum >= start
        snapNum = snapNum[ii]
        mstar = mstar[ii]
        fgas = fgas[ii]
        snapNumNext = snapNumNext[ii]
        mstarNext = mstarNext[ii]
        fgasNext = fgasNext[ii]
        R11 = 0.3
        R9 = 0.2
        R = normRadius(mstar, R11=R11, R9=R9)
        Rnext = normRadius(mstarNext, R11=R11, R9=R9)
        norm = colors.Normalize(vmin=0.0, vmax=0.2)
        color_f = norm(fgas)
        color_n = norm(fgasNext)
        cmap = plt.get_cmap('seismic_r')
        plt.clf()
        fig = plt.figure(figsize=[ii.sum()/3.-4.0, 5./3.])
        fig.subplots_adjust(left=0.01, bottom=0.15, right=0.99,
                            top=0.95, hspace=0.2, wspace=0.4)
        ax = fig.add_subplot(111)
        ax.set_aspect(1)
        ax.set_xlim(snapNum[-1]-2, snapNum[0]+2)
        ax.set_ylim(0, 5)
        util_fig.set_labels(ax)
        for i in range(ii.sum()):
            circle_f = mpatches.Circle((snapNum[i], 1), R[i],
                                       facecolor=cmap(color_f[i]),
                                       edgecolor=None, alpha=1, zorder=1)
            ax.add_patch(circle_f)
            if i > 0:
                ax.annotate('', xy=(snapNum[i-1], 1), xycoords='data',
                            xytext=(snapNum[i], 1), textcoords='data',
                            arrowprops=dict(arrowstyle='-'), zorder=0)
            if not np.isnan(Rnext[i]):
                circle_n = mpatches.Circle((snapNumNext[i], 2.5), Rnext[i],
                                           facecolor=cmap(color_n[i]),
                                           edgecolor=None,
                                           alpha=1, zorder=1)
                ax.add_patch(circle_n)
                ax.annotate('', xy=(snapNum[i]+1, 1), xycoords='data',
                            xytext=(snapNumNext[i], 2.5), textcoords='data',
                            arrowprops=dict(arrowstyle='-'), zorder=0)

        circle = mpatches.Circle((snapNum[-1], 3.5),
                                 normRadius(12.0, R11=R11, R9=R9),
                                 facecolor='k', edgecolor=None, alpha=1,
                                 zorder=1)
        ax.text(snapNum[-1]+0.6, 4.3, '12', fontproperties=util_fig.ticks_font)
        ax.add_patch(circle)

        circle = mpatches.Circle((snapNum[-3], 3.5),
                                 normRadius(11.0, R11=R11, R9=R9),
                                 facecolor='k', edgecolor=None, alpha=1,
                                 zorder=1)
        ax.text(snapNum[-3]+0.6, 4.3, '11', fontproperties=util_fig.ticks_font)
        ax.add_patch(circle)

        circle = mpatches.Circle((snapNum[-5], 3.5),
                                 normRadius(10.0, R11=R11, R9=R9),
                                 facecolor='k', edgecolor=None, alpha=1,
                                 zorder=1)
        ax.text(snapNum[-5]+0.6, 4.3, '10', fontproperties=util_fig.ticks_font)
        ax.add_patch(circle)
        circle = mpatches.Circle((snapNum[-7], 3.5),
                                 normRadius(9.0, R11=R11, R9=R9),
                                 facecolor='k', edgecolor=None, alpha=1,
                                 zorder=1)
        ax.text(snapNum[-7]+0.6, 4.3, '9', fontproperties=util_fig.ticks_font)
        ax.add_patch(circle)
        ax.set_yticklabels([])
        ax.axvline(85.0, linestyle='dashed', color='r', lw=3.0)
        ax.axvline(103.0, linestyle='dashed', color='r', lw=3.0)
        axc = fig.add_axes([0.3, 0.79, 0.35, 0.11])
        util_fig.set_labels(axc)
        axc.set_yticklabels([])
        colorbar.ColorbarBase(axc, cmap=cmap, norm=norm,
                              orientation='horizontal')
        fig.savefig('{}/mpbTree.pdf'.format(outpath), dpi=500)
