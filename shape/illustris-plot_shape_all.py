#!/usr/bin/env python
import numpy as np
import matplotlib
matplotlib.use('Agg')
from JAM.utils import util_fig
import matplotlib.pyplot as plt
from matplotlib import colors
from optparse import OptionParser
import sys
import pickle
import util_illustris as ui
from illustris_shape_util import vect2angle
from illustris_shape_util import plot_shape_profile
from illustris_shape_util import plot_shape_dot
from illustris_shape_util import plot_orientation_profile
from illustris_shape_util import plot_orientation_dot
util_fig.label_font.set_size(15)
util_fig.ticks_font.set_size(10)
util_fig.text_font.set_size(10)
sauron = util_fig.sauron


def plot_all(Rstar, axisRatiosStar, hmr_star, eigenVectorsStar,
             Rdark, axisRatiosDark, hmr_dark, eigenVectorsDark,
             path='.'):
    # plot shape profile for star
    fig = plt.figure(figsize=(4, 3))
    fig.subplots_adjust(left=0.155, bottom=0.15, right=0.97, top=0.97)
    ax = fig.add_subplot(111)
    plot_shape_profile(Rstar, axisRatiosStar, axes=ax, color='r')
    ax.axvline(hmr_star, color='c', lw=2.0)
    ax.axvline(2.0*hmr_star, color='c', lw=2.0)
    ax.set_xlim([0.0, 3*hmr_star])
    ax.set_xlabel('R [kpc]', fontproperties=util_fig.label_font)
    ax.set_ylabel('Axis ratios', fontproperties=util_fig.label_font)
    util_fig.set_labels(ax)
    fig.savefig('{}/shape/shape_star.png'.format(path), dpi=500)

    # plot shape dot for star
    fig = plt.figure(figsize=(4, 4))
    fig.subplots_adjust(left=0.14, bottom=0.12, right=0.97, top=0.97)
    ax = fig.add_subplot(111)
    plot_shape_dot(Rstar, axisRatiosStar, axes=ax, color='r')
    # plot_shape_dot(Rdark, axisRatiosDark, axes=ax, color='k')
    ax.set_xlabel('b/a', fontproperties=util_fig.label_font)
    ax.set_ylabel('c/a', fontproperties=util_fig.label_font)
    lim = np.array([0, 1])
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    length = 0.2
    x = np.linspace(0.0, 1.0, 10)
    y = x - length
    ax.plot(lim, lim, '--k', lw=2.5)
    ax.plot(x, y, '--k', lw=2.5)
    ax.plot([1-length, 1-length], [0.0, 1-length], '--k', lw=2.5, zorder=1)
    util_fig.set_labels(ax)
    fig.savefig('{}/shape/shape_star_dot.png'.format(path), dpi=500)


def plot_all_main(path):
    flist = np.genfromtxt('{}/select.list'.format(path), dtype=[('snap', 'S15'), ('ID', 'S15')])
    fig_dot = plt.figure(figsize=(4, 4))
    fig_dot.subplots_adjust(left=0.125, bottom=0.105, right=0.97, top=0.97)
    ax_dot = fig_dot.add_subplot(111)
    index = 0.3
    norm = colors.Normalize(vmin=(0.0), vmax=(3.0**index))
    text_x, text_y = 0.1, 0.9
    for i in range(len(flist)):
        fname = '{}-{}'.format(flist['snap'][i], flist['ID'][i])
        z = ui.snap2z(int(flist['snap'][i][4:]))
        snap = flist['snap'][i][4:]
        try:
            with open('{}/{}/shape/shape.dat'.format(path, fname), 'rb') as f:
                data = pickle.load(f)
            Rstar = data['Rstar']
            hmr_star = data['hmr_star']
            axisRatiosStar = data['axisRatiosStar']
            # eigenVectorsStar = data['eigenVectorsStar']
            # Rdark = data['Rdark']
            # hmr_dark = data['hmr_dark']
            # axisRatiosDark = data['axisRatiosDark']
            # eigenVectorsDark = data['eigenVectorsDark']
            plot_shape_dot(Rstar, axisRatiosStar, axes=ax_dot, color=sauron(norm(z**index)))
            n = 4
            shift_x = 0.1 * (i%n)
            shift_y = 0.1 * (i//n)
            ax_dot.text(text_x+shift_x, text_y-shift_y, snap, color=sauron(norm(z**index)),
                        transform=ax_dot.transAxes, fontproperties=util_fig.text_font)
        except:
            pass
    ax_dot.set_xlabel('b/a', fontproperties=util_fig.label_font)
    ax_dot.set_ylabel('c/a', fontproperties=util_fig.label_font)
    lim_dot = np.array([0, 1])
    ax_dot.set_xlim(lim_dot)
    ax_dot.set_ylim(lim_dot)
    length = 0.2
    x = np.linspace(0.0, 1.0, 10)
    y = x - length
    ax_dot.plot(x, y, '--k', lw=1.5, zorder=0)
    ax_dot.plot([1-length, 1-length], [0.0, 1-length], '--k', lw=1.5, zorder=0)
    ax_dot.plot(lim_dot, lim_dot, '--k', lw=1.5, zorder=0)
    util_fig.set_labels(ax_dot)
    fig_dot.savefig('{}/shape_star_dot.pdf'.format(path), dpi=500)
    '''
    with open('{}/shape/shape.dat'.format(path), 'rb') as f:
        data = pickle.load(f)
    Rstar = data['Rstar']
    hmr_star = data['hmr_star']
    axisRatiosStar = data['axisRatiosStar']
    eigenVectorsStar = data['eigenVectorsStar']
    Rdark = data['Rdark']
    hmr_dark = data['hmr_dark']
    axisRatiosDark = data['axisRatiosDark']
    eigenVectorsDark = data['eigenVectorsDark']
    '''

if __name__ == '__main__':
    parser = OptionParser()
    (options, args) = parser.parse_args()
    if len(args) != 1:
        print 'Error - please provide a folder name'
        sys.exit(1)
    plot_all_main(args[0])
