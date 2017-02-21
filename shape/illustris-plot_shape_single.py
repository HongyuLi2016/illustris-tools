#!/usr/bin/env python
import numpy as np
import matplotlib
matplotlib.use('Agg')
from JAM.utils import util_fig
import matplotlib.pyplot as plt
from optparse import OptionParser
import sys
import pickle
from illustris_shape_util import vect2angle
from illustris_shape_util import plot_shape_profile
from illustris_shape_util import plot_shape_dot
from illustris_shape_util import plot_orientation_profile
from illustris_shape_util import plot_orientation_dot
util_fig.label_font.set_size(15)
util_fig.ticks_font.set_size(10)


def plot_single(Rstar, axisRatiosStar, hmr_star, eigenVectorsStar,
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
    ax.plot(lim, lim, '--k', lw=1)
    util_fig.set_labels(ax)
    fig.savefig('{}/shape/shape_star_dot.png'.format(path), dpi=500)

    # plot shape profile for dark + star
    fig = plt.figure(figsize=(4, 3))
    fig.subplots_adjust(left=0.155, bottom=0.15, right=0.97, top=0.97)
    ax = fig.add_subplot(111)
    plot_shape_profile(Rstar, axisRatiosStar, axes=ax, color='r')
    plot_shape_profile(Rdark, axisRatiosDark, axes=ax, color='k')
    ax.axvline((hmr_dark), color='c', lw=2.0)
    ax.axvline(2.0*hmr_dark, color='c', lw=2.0)
    ax.set_xlim([0.0, 2.5*hmr_dark])
    ax.set_xlabel('R [kpc]', fontproperties=util_fig.label_font)
    ax.set_ylabel('Axis ratios', fontproperties=util_fig.label_font)
    util_fig.set_labels(ax)
    fig.savefig('{}/shape/shape_dark.png'.format(path), dpi=500)

    # plot shape dot for star + dark
    fig = plt.figure(figsize=(4, 4))
    fig.subplots_adjust(left=0.14, bottom=0.12, right=0.97, top=0.97)
    ax = fig.add_subplot(111)
    plot_shape_dot(Rstar, axisRatiosStar, axes=ax, color='r')
    plot_shape_dot(Rdark, axisRatiosDark, axes=ax, color='k')
    ax.set_xlabel('b/a', fontproperties=util_fig.label_font)
    ax.set_ylabel('c/a', fontproperties=util_fig.label_font)
    lim = np.array([0, 1])
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.plot(lim, lim, '--k', lw=1)
    util_fig.set_labels(ax)
    fig.savefig('{}/shape/shape_dark_dot.png'.format(path), dpi=500)

    # plot orientation profile for star
    fig, axes = plt.subplots(2, 3)
    plot_orientation_profile(Rstar, eigenVectorsStar,
                             axes=axes, color='r')
    axes[1, 1].set_xlabel('R [kpc]', fontproperties=util_fig.label_font)
    axes[0, 0].set_ylabel('theta', fontproperties=util_fig.label_font)
    axes[1, 0].set_ylabel('phi', fontproperties=util_fig.label_font)
    for ax in axes.reshape(-1):
        util_fig.set_labels(ax)
    for ax in axes[0, :]:
        ax.set_ylim([0, 90])
    for ax in axes[1, :]:
        ax.set_ylim([-180, 180])
    fig.savefig('{}/shape/orientation_star.png'.format(path), dpi=500)

    # plot orientation profile for star + dark
    fig, axes = plt.subplots(2, 3)
    plot_orientation_profile(Rstar, eigenVectorsStar,
                             axes=axes, color='r')
    plot_orientation_profile(Rdark, eigenVectorsDark,
                             axes=axes, color='k')
    axes[1, 1].set_xlabel('R [kpc]', fontproperties=util_fig.label_font)
    axes[0, 0].set_ylabel('theta', fontproperties=util_fig.label_font)
    axes[1, 0].set_ylabel('phi', fontproperties=util_fig.label_font)
    for ax in axes.reshape(-1):
        util_fig.set_labels(ax, xrotate=45)
    for ax in axes[0, :]:
        ax.set_ylim([0, 90])
    for ax in axes[1, :]:
        ax.set_ylim([-180, 180])
    fig.savefig('{}/shape/orientation_dark.png'.format(path), dpi=500)

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    plot_orientation_dot(Rstar, eigenVectorsStar,
                         axes=axes, color='r')
    # plot_orientation_dot(Rdark, eigenVectorsDark, fig=fig,
    #                      axes=axes, color='k')
    axes[1].set_xlabel('phi', fontproperties=util_fig.label_font)
    axes[0].set_ylabel('theta', fontproperties=util_fig.label_font)
    for ax in axes[:]:
        ax.set_ylim([0, 90])
        ax.set_xlim([-180, 180])
        util_fig.set_labels(ax)
    fig.subplots_adjust(left=0.04, bottom=0.13, right=0.97, top=0.97)
    fig.savefig('{}/shape/orientation_star_dot.png'.format(path), dpi=500)

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    plot_orientation_dot(Rstar, eigenVectorsStar,
                         axes=axes, color='r')
    plot_orientation_dot(Rdark, eigenVectorsDark,
                         axes=axes, color='k')
    axes[1].set_xlabel('phi', fontproperties=util_fig.label_font)
    axes[0].set_ylabel('theta', fontproperties=util_fig.label_font)
    for ax in axes[:]:
        ax.set_ylim([0, 90])
        ax.set_xlim([-180, 180])
        util_fig.set_labels(ax)
    fig.subplots_adjust(left=0.04, bottom=0.13, right=0.97, top=0.97)
    fig.savefig('{}/shape/orientation_dark_dot.png'.format(path), dpi=500)


def plot_all():
    pass


def plot_single_main(path):
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
    plot_single(Rstar, axisRatiosStar, hmr_star, eigenVectorsStar,
                Rdark, axisRatiosDark, hmr_dark, eigenVectorsDark,
                path=path)


if __name__ == '__main__':
    parser = OptionParser()
    (options, args) = parser.parse_args()
    if len(args) != 1:
        print 'Error - please provide a folder name'
        sys.exit(1)
    plot_single_main(args[0])
