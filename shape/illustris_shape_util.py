#!/usr/bin/env python
import numpy as np
from JAM.utils import util_fig
util_fig.label_font.set_size(15)
util_fig.ticks_font.set_size(10)


def vect2angle(eigenVector):
    '''
    convert eigen vector to angle theta, phi
    '''
    if eigenVector[2] < 0.0:
        eigenVector *= -1.0
    eigenVector /= np.sum(eigenVector**2)**0.5  # normalize
    theta = np.arccos(eigenVector[2])
    if eigenVector[2] == 1.0:
        phi = 0.0
    else:
        cosphi = eigenVector[0] / np.sqrt(1 - eigenVector[2]**2)
        phi_tem = np.arccos(cosphi)
        if eigenVector[1] >= 0:
            phi = phi_tem
        else:
            phi = -phi_tem
    return theta, phi


def plot_shape_profile(R, axisRatios, axes=None, color='k'):
    axes.plot(R, axisRatios[:, 0], '-', lw=3, color=color)
    axes.plot(R, axisRatios[:, 1], '--', lw=3, color=color)


def plot_shape_dot(R, axisRatios, axes=None, color='k',
                   markersize=2.0):
    size = (R - R.min())/(R.max() - R.min()) * 3 + markersize
    for i in range(len(R)):
        axes.plot(axisRatios[i, 0], axisRatios[i, 1], 'o',
                  markerfacecolor=color, markeredgecolor='none', lw=0.05,
                  markersize=size[i])


def plot_orientation_profile(R, eigenVectors, axes=None, color='k'):
    theta = np.zeros([len(R), 3])
    phi = np.zeros([len(R), 3])
    for i in range(len(R)):
        theta[i, 0], phi[i, 0] = vect2angle(eigenVectors[i, 0, :])
        theta[i, 1], phi[i, 1] = vect2angle(eigenVectors[i, 1, :])
        theta[i, 2], phi[i, 2] = vect2angle(eigenVectors[i, 2, :])
    axes[0, 0].plot(R, np.degrees(theta[:, 0]), color=color, lw=2)
    axes[0, 1].plot(R, np.degrees(theta[:, 1]), color=color, lw=2)
    axes[0, 2].plot(R, np.degrees(theta[:, 2]), color=color, lw=2)
    axes[1, 0].plot(R, np.degrees(phi[:, 0]), color=color, lw=2)
    axes[1, 1].plot(R, np.degrees(phi[:, 1]), color=color, lw=2)
    axes[1, 2].plot(R, np.degrees(phi[:, 2]), color=color, lw=2)


def plot_orientation_dot(R, eigenVectors, axes=None, color='k',
                         markersize=2.0, setlim=True):
    theta = np.zeros([len(R), 3])
    phi = np.zeros([len(R), 3])
    size = (R - R.min())/(R.max() - R.min()) * 3 + markersize
    for i in range(len(R)):
        theta[i, 0], phi[i, 0] = vect2angle(eigenVectors[i, 0, :])
        theta[i, 1], phi[i, 1] = vect2angle(eigenVectors[i, 1, :])
        theta[i, 2], phi[i, 2] = vect2angle(eigenVectors[i, 2, :])
        axes[0].plot(np.degrees(phi[i, 0]), np.degrees(theta[i, 0]),
                     'o', markerfacecolor=color, markeredgecolor='none',
                     markersize=size[i])
        axes[1].plot(np.degrees(phi[i, 1]), np.degrees(theta[i, 1]),
                     'o', markerfacecolor=color, markeredgecolor='none',
                     markersize=size[i])
        axes[2].plot(np.degrees(phi[i, 2]), np.degrees(theta[i, 2]),
                     'o', markerfacecolor=color, markeredgecolor='none',
                     markersize=size[i])
    if setlim:
        axes[0].set_xlim([-180, 180])
        axes[0].set_ylim([0, 90])
        axes[1].set_xlim([-180, 180])
        axes[1].set_ylim([0, 90])
        axes[2].set_xlim([-180, 180])
        axes[2].set_ylim([0, 90])
