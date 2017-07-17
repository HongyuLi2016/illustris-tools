import numpy as np
import h5py
import numpy.linalg as LA
import requests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
from JAM.utils import util_fig
from scipy.stats import binned_statistic_2d
import re
import warnings

warnings.simplefilter("ignore")
util_fig.ticks_font.set_size('x-small')

# global parameters
boxsize_img = 80.0
boxsize_dark = 150
scale_img = 0.5  # pixel2kpc
scale_dark = 1.0
boxsize_ifu = 80.0
scale_ifu = 1.00
kpc2arcsec = 1.612
h0 = 0.704
massUnit = 1e10/h0
MdarkPart = 6.3e6


def snap2z(snapNum):
    snapNum = np.atleast_1d(snapNum)
    z = [4.67730473e+01, 4.45622038e+01, 4.24536738e+01, 4.06395569e+01,
         3.87125594e+01, 3.68747395e+01, 3.51219704e+01, 3.36139397e+01,
         3.20120740e+01, 3.04843396e+01, 2.90273057e+01, 2.76377005e+01,
         2.64421253e+01, 2.51721572e+01, 2.39609608e+01, 2.28058162e+01,
         2.18119639e+01, 2.07562707e+01, 1.97494329e+01, 1.87891896e+01,
         1.79630246e+01, 1.70854528e+01, 1.62484933e+01, 1.54502666e+01,
         1.47634960e+01, 1.40339921e+01, 1.33382483e+01, 1.26747021e+01,
         1.20418635e+01, 1.14973880e+01, 1.09190332e+01, 1.03674436e+01,
         9.99659047e+00, 9.84138044e+00, 9.38877127e+00, 9.00233985e+00,
         8.90799919e+00, 8.44947629e+00, 8.01217295e+00, 7.59510715e+00,
         7.23627607e+00, 7.00541705e+00, 6.85511726e+00, 6.49159775e+00,
         6.14490120e+00, 6.01075740e+00, 5.84661375e+00, 5.52976581e+00,
         5.22758097e+00, 4.99593347e+00, 4.93938066e+00, 4.66451770e+00,
         4.42803374e+00, 4.00794511e+00, 3.70877426e+00, 3.49086137e+00,
         3.28303306e+00, 3.08482264e+00, 3.00813107e+00, 2.89578501e+00,
         2.73314262e+00, 2.57729027e+00, 2.44422570e+00, 2.31611074e+00,
         2.20792547e+00, 2.10326965e+00, 2.00202814e+00, 1.90408954e+00,
         1.82268925e+00, 1.74357057e+00, 1.66666956e+00, 1.60423452e+00,
         1.53123903e+00, 1.47197485e+00, 1.41409822e+00, 1.35757667e+00,
         1.30237846e+00, 1.24847261e+00, 1.20625808e+00, 1.15460271e+00,
         1.11415056e+00, 1.07445789e+00, 1.03551045e+00, 9.97294226e-01,
         9.87852811e-01, 9.50531352e-01, 9.23000816e-01, 8.86896938e-01,
         8.51470901e-01, 8.16709979e-01, 7.91068249e-01, 7.57441373e-01,
         7.32636182e-01, 7.00106354e-01, 6.76110411e-01, 6.44641841e-01,
         6.21428745e-01, 5.98543288e-01, 5.75980845e-01, 5.46392183e-01,
         5.24565820e-01, 5.03047523e-01, 4.81832943e-01, 4.60917794e-01,
         4.40297849e-01, 4.19968942e-01, 3.99926965e-01, 3.80167867e-01,
         3.60687657e-01, 3.47853842e-01, 3.28829724e-01, 3.10074120e-01,
         2.91583240e-01, 2.73353347e-01, 2.61343256e-01, 2.43540182e-01,
         2.25988386e-01, 2.14425036e-01, 1.97284182e-01, 1.80385262e-01,
         1.69252033e-01, 1.52748769e-01, 1.41876204e-01, 1.25759332e-01,
         1.09869940e-01, 9.94018026e-02, 8.38844308e-02, 7.36613847e-02,
         5.85073228e-02, 4.85236300e-02, 3.37243719e-02, 2.39744284e-02,
         9.52166697e-03, 2.22044605e-16]
    number = [0,   1,   2,   3,   4,   5,   6,   7,   8,   9,
              10,  11,  12,  13,  14,  15,  16,  17,  18,  19,
              20,  21,  22,  23,  24,  25,  26,  27,  28,  29,
              30,  31,  32,  33,  34,  35,  36,  37,  38,  39,
              40,  41,  42,  43,  44,  45,  46,  47,  48,  49,
              50,  51,  52,  54,  56,  57,  58,  59,  60,  61,
              62,  63,  64,  65,  66,  67,  68,  69,  70,  71,
              72,  73,  74,  75,  76,  77,  78,  79,  80,  81,
              82,  83,  84,  85,  86,  87,  88,  89,  90,  91,
              92,  93,  94,  95,  96,  97,  98,  99,  100, 101,
              102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
              112, 113, 114, 115, 116, 117, 118, 119, 120, 121,
              122, 123, 124, 125, 126, 127, 128, 129, 130, 131,
              132, 133, 134, 135]
    z = np.asarray(z)
    number = np.asarray(number)
    ii = number == snapNum
    if ii.sum() == 1:
        return z[ii][0]
    else:
        return np.nan


def extract_id(fname):
    try:
        pat = re.compile(r'snap([-+]?\d+)-subhalo([-+]?\d+)')
        search = pat.search(fname)
        snapNum = search.group(1)
        subhaloID = search.group(2)
    except:
        pat = re.compile(r'snap([-+]?\d+)-halo([-+]?\d+)')
        search = pat.search(fname)
        snapNum = search.group(1)
        subhaloID = search.group(2)
    return int(snapNum), int(subhaloID)


def vel_center(xpart, vpart, mpart=None, R=15.0):
    if mpart is None:
        mpart = np.ones(xpart.shape[0], dtype=float)
    r = np.sqrt(np.sum(xpart**2, axis=1))
    iIn = r < R
    return np.average(vpart[iIn], weights=mpart[iIn], axis=0)


def find_center(xpart, mpart=None, percent=20.0, imax=100):
    if mpart is None:
        mpart = np.ones(xpart.shape[0], dtype=float)
    xcenter_old = np.average(xpart, axis=0, weights=mpart)
    deltaD = 10.0
    # iteratively find out the mass center
    i = 0
    while deltaD > 0.01:
        r = np.sqrt((xpart[:, 0]-xcenter_old[0])**2 +
                    (xpart[:, 1]-xcenter_old[1])**2 +
                    (xpart[:, 2]-xcenter_old[2])**2)
        iIn = r < np.percentile(r, percent)
        xcenter_new = \
            np.average(xpart[iIn, :], axis=0, weights=mpart[iIn])
        deltaD = np.sum((xcenter_new - xcenter_old)**2)**0.5
        xcenter_old = xcenter_new
        i += 1
        if i > imax:
            print('Warning - find center reach maximum iteration')
            break
    return xcenter_new


def mock_img(x1, x2, L=None, bins=100, Xrange=None):
    if L is None:
        L = np.ones(len(x1))
    per = [5, 95]
    if Xrange is None:
        range1 = np.percentile(x1, per)
        range2 = np.percentile(x2, per)
        Xrange = np.array([range1, range2])
    H, x1edge, x2edge = np.histogram2d(x1, x2, bins=bins, range=Xrange,
                                       weights=L)
    return H, x1edge, x2edge


def mock_ifu(x1, x2, vel, bins=100, Xrange=None):
    per = [5, 95]
    if Xrange is None:
        range1 = np.percentile(x1, per)
        range2 = np.percentile(x2, per)
        Xrange = np.array([range1, range2])
    vmap, x1edge, x2edge, binNum = \
        binned_statistic_2d(x1, x2, vel, bins=bins,
                            statistic='mean', range=Xrange)
    return vmap, x1edge, x2edge


def read_cutout(fname, PartType=4, key='Coordinates', z=0.0):
    '''
    keys:
      dm - Coordinates, Velocities
      stars - Coordinates, Velocities, GFM_Metallicity,
              GFM_StellarFormationTime, GFM_StellarPhotometrics, Masses,
    units are convert to physical units (i.e. no h0 and scale factor a)
    '''
    scale_factor = 1.0 / (1+z)
    with h5py.File(fname, 'r') as f:
        Ptype = 'PartType{}'.format(PartType)
        data = f[Ptype][key][:]
    if key == 'Coordinates':
        # distance in kpc
        data *= scale_factor/h0
    elif key == 'Velocities':
        # velocity in km/s
        data *= np.sqrt(scale_factor)
    elif key == 'Masses':
        # mass in M_sun
        data /= h0
    return data


class cutout(object):
    '''
    read cutout info, all the values are in original units
    '''
    def __init__(self, fname):
        with h5py.File(fname, 'r') as f:
            data_dark = f['PartType1']
            data_star = f['PartType4']
            self.numDark = len(data_dark['Coordinates'][:][:, 0])
            self.numStar = len(data_star['Coordinates'][:][:, 0])
            self.Mdark = np.log10(self.numDark * MdarkPart)
            self.Mstar = np.log10(data_star['Masses'][:].sum()*massUnit)
            self.xpart_star = data_star['Coordinates'][:]
            self.vpart_star = data_star['Velocities'][:]
            self.xpart_dark = data_dark['Coordinates'][:]
            self.vpart_dark = data_dark['Velocities'][:]
            self.mass_star = data_star['Masses'][:] * massUnit
            self.mass_dark = MdarkPart * \
                np.ones_like(data_dark['Coordinates'][:][:, 0], dtype=float)


def read_Dandan(snap, key, path='/share/data/D/lhy/paper4/Dandan_catalogue'):
    if snap == 135:
        if key == 'subfindID':
            with h5py.File('{}/snap135.hdf5'.format(path), 'r') as f:
                return (f['FID'][:] + f['SID'][:]).astype(int)[:, 0]
        else:
            with h5py.File('{}/snap135.hdf5'.format(path), 'r') as f:
                return f[key][:][:, 0]
    else:
        with h5py.File('{}/snaps.hdf5'.format(path), 'r') as f:
            snapID = 'Snapshot_{}'.format(snap)
            return f[snapID][key][:][:, 0]


def cutout_vel_los(xpart, vpart, mpart, bins=100, box=None, magrange=5,
                   xmark=None, **kwargs):
    '''
    plot image and los vel
    '''
    projection = [(0, 1, 2, 0), (0, 2, 1, 1), (1, 2, 0, 2)]
    axis = ['x', 'y', 'z']
    fig, axes = plt.subplots(2, 3, figsize=(9, 6))
    fig.subplots_adjust(left=0.12, bottom=0.1, right=0.9,
                        top=0.95, hspace=0.2, wspace=0.4)
    if box is None:
        box = np.percentile(xpart, [5, 95], axis=0)
    for i, j, los, k in projection:
        Xrange = [[box[0, i], box[1, i]],
                  [box[0, j], box[1, j]]]
        # make mock images
        img, x1edge, x2edge = \
            mock_img(xpart[:, i], xpart[:, j], L=mpart,
                     bins=bins, Xrange=Xrange)
        extent = np.array([x1edge[0], x1edge[-1],
                           x2edge[0], x2edge[-1]])
        vel, x1edge, x2edge = \
            mock_ifu(xpart[:, i], xpart[:, j], vpart[:, los], bins=100,
                     Xrange=Xrange)
        axes[0, k].imshow(np.rot90(np.log10(img)), extent=extent,
                          cmap='gray')
        # over plot contours
        x1 = np.linspace(Xrange[0][0], Xrange[0][1], bins)
        x2 = np.linspace(Xrange[1][0], Xrange[1][1], bins)
        X1, X2 = np.meshgrid(x1, x2, indexing='ij')
        good = ~np.isnan(img)
        Imax = np.percentile(img[good].max(), 99.5)
        levels = Imax * 10**(-0.4*np.arange(0, magrange, 0.5))
        levels = levels[::-1]
        axes[0, k].contour(X1, X2, img, levels=levels,
                           linestyles='solid', **kwargs)
        # mark some points
        if xmark is not None:
            xmark = np.atleast_2d(xmark)
            axes[0, k].plot(xmark[:, i], xmark[:, j], '+r', markersize=5,)

        good = ~np.isnan(vel)
        vmax = np.percentile(abs(vel[good]), 95)
        norm = colors.Normalize(vmin=-vmax, vmax=vmax)
        axes[1, k].imshow(np.rot90(vel), extent=extent,
                          cmap=util_fig.sauron, norm=norm)
        # add text
        axes[0, k].text(0.1, 0.9, '{}{}'.format(axis[i], axis[j]),
                        transform=axes[0, k].transAxes,
                        fontproperties=util_fig.label_font)
        util_fig.add_colorbar(axes[1, k], util_fig.sauron, norm,
                              fig=fig, width=0.01)
        util_fig.set_labels(axes[0, k], xrotate=45)
        util_fig.set_labels(axes[1, k], xrotate=45)
        axes[1, k].set_xlabel(r'$\rm cKpc/h_0$',
                              fontproperties=util_fig.label_font)
    axes[0, 0].set_ylabel(r'$\rm cKpc/h_0$',
                          fontproperties=util_fig.label_font)
    axes[1, 0].set_ylabel(r'$\rm cKpc/h_0$',
                          fontproperties=util_fig.label_font)
    return fig, axes


def read_center_mark(fname):
    cmark = []
    with open(fname) as f:
        num_mark = int(f.readline().splitlines()[0])
        tem_cmark = f.readline().splitlines()[0].split()
        for i in range(num_mark):
            cmark.append(tem_cmark[i])
        pos = np.genfromtxt(fname, skip_header=2)
    return pos, cmark


def cutout_vel_vector(xpart, vpart, mpart, bins=100, ifu_bins=50, box=None,
                      magrange=5, linewidths=0.3, contoursColor='c',
                      xmark=None, xmarkColor=None, **kwargs):
    '''
    plot img and vector v feild
    box: 2*3 array, boundary of the box within witch img and vmap will be
         calculated
    '''
    projection = [(0, 1, 2, 0), (0, 2, 1, 1), (1, 2, 0, 2)]
    axis = ['x', 'y', 'z']

    fig, axes = plt.subplots(2, 3, figsize=(9, 6))
    fig.subplots_adjust(left=0.12, bottom=0.1, right=0.9,
                        top=0.95, hspace=0.2, wspace=0.4)
    if box is None:
        box = np.percentile(xpart, [5, 95], axis=0)
    for i, j, los, k in projection:
        Xrange = [[box[0, i], box[1, i]],
                  [box[0, j], box[1, j]]]
        # make mock images
        img, x1edge, x2edge = \
            mock_img(xpart[:, i], xpart[:, j], L=mpart,
                     bins=bins, Xrange=Xrange)
        extent = np.array([x1edge[0], x1edge[-1],
                           x2edge[0], x2edge[-1]])

        axes[0, k].imshow(np.rot90(np.log10(img)), extent=extent,
                          cmap='gray')
        axes[1, k].imshow(np.rot90(np.log10(img)), extent=extent,
                          cmap='gray', zorder=0)

        # over plot contours
        x1 = np.linspace(Xrange[0][0], Xrange[0][1], bins)
        x2 = np.linspace(Xrange[1][0], Xrange[1][1], bins)
        X1, X2 = np.meshgrid(x1, x2, indexing='ij')
        good = ~np.isnan(img)
        Imax = np.percentile(img[good].max(), 99.5)
        levels = Imax * 10**(-0.4*np.arange(0, magrange, 0.5))
        levels = levels[::-1]
        axes[0, k].contour(X1, X2, img, levels=levels,
                           linestyles='solid', linewidths=linewidths,
                           colors=contoursColor)

        # mark some points
        if xmark is not None:
            xmark = np.atleast_2d(xmark)
            if xmarkColor is None:
                xmarkColor = ['r'] * xmark.shape[0]
            for cnum in range(xmark.shape[0]):
                axes[0, k].plot(xmark[cnum, i], xmark[cnum, j], '+',
                                color=xmarkColor[cnum], markersize=5,)
        # over plot velocity feild
        vel_x1, v1edge_star, v2edge_star = \
            mock_ifu(xpart[:, i], xpart[:, j], vpart[:, i], bins=ifu_bins,
                     Xrange=Xrange)
        vel_x2, v1edge_star, v2edge_star = \
            mock_ifu(xpart[:, i], xpart[:, j], vpart[:, j], bins=ifu_bins,
                     Xrange=Xrange)
        vel_los, x1edge_star, x2edge_star = \
            mock_ifu(xpart[:, i], xpart[:, j], vpart[:, los], bins=ifu_bins,
                     Xrange=Xrange)

        # remove the system velocity
        good = ~np.isnan(vel_los)
        vmax = np.percentile(abs(vel_los[good]), 95)
        norm = colors.Normalize(vmin=-vmax, vmax=vmax)

        x1 = np.linspace(Xrange[0][0], Xrange[0][1], ifu_bins)
        x2 = np.linspace(Xrange[1][0], Xrange[1][1], ifu_bins)
        X1, X2 = np.meshgrid(x1, x2, indexing='ij')

        color = vel_los
        axes[1, k].quiver(X1, X2, vel_x1, vel_x2, color, norm=norm,
                          cmap=util_fig.sauron, units='inches', width=0.01,
                          pivot='tip', zorder=1, **kwargs)

        # add text
        axes[1, k].text(0.1, 0.9, '{}{}'.format(axis[i], axis[j]),
                        transform=axes[0, k].transAxes,
                        fontproperties=util_fig.label_font)
        util_fig.set_labels(axes[0, k], xrotate=45)
        util_fig.set_labels(axes[1, k], xrotate=45)
        axes[1, k].set_xlabel(r'$\rm cKpc/h_0$',
                              fontproperties=util_fig.label_font)
    axes[0, 0].set_ylabel(r'$\rm cKpc/h_0$',
                          fontproperties=util_fig.label_font)
    axes[1, 0].set_ylabel(r'$\rm cKpc/h_0$',
                          fontproperties=util_fig.label_font)
    return fig, axes


def get_spin(x, v, mpart=None, Rb=20):
    r = np.sqrt(x[:, 0]**2 + x[:, 1]**2 + x[:, 2]**2)
    good = r < Rb
    if mpart is None:
        mpart = np.ones(len(good))
    mpart = mpart[good]
    x = x[good, :]
    v = v[good, :]
    Spin = np.dot(mpart.reshape(1, -1), np.cross(x, v, axis=1))[0]
    spinValue = np.linalg.norm(Spin)
    spinVector = Spin / spinValue
    Mpart = mpart.sum()
    return spinValue, Mpart, spinVector


def get_shape(x, mpart, Rb=20., decrease=True):
    # Rb=20kpc, within which the shape is calcualted
    s = 1.
    q = 1.

    Tiv = np.array([[1, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1]])
    # tmp = Tiv

    order = [2, 1, 0]

    Vei = np.zeros((3, 3))

    dq = 10000.
    ds = 10000.

    while (dq > 0.01 or ds > 0.01):
        # in eigenvector coordinates
        y = np.transpose(np.dot(Tiv, np.transpose(x)))
        rn0 = np.sqrt(np.power(y[:, order[2]], 2.) +
                      np.power(y[:, order[1]], 2.)/q/q +
                      np.power(y[:, order[0]], 2.)/s/s)
        ind = np.where(rn0 < Rb)[0]
        # Np = ind.shape[0]

        y1 = y[ind, 0]
        y2 = y[ind, 1]
        y3 = y[ind, 2]
        rn = rn0[ind]

        I11 = np.sum(y1*y1/np.power(rn, 2))
        I22 = np.sum(y2*y2/np.power(rn, 2))
        I33 = np.sum(y3*y3/np.power(rn, 2))
        I12 = np.sum(y1*y2/np.power(rn, 2))
        I13 = np.sum(y1*y3/np.power(rn, 2))
        I23 = np.sum(y2*y3/np.power(rn, 2))

        II = [[I11, I12, I13],
              [I12, I22, I23],
              [I13, I23, I33]]

        D, A = LA.eig(II)
        # print 'eigenvalues'
        # print D
        # print  'eigenvectors'
        # print A
        order = np.argsort(D)  # a=order2,b=order1,c=order0
        la = np.sqrt(D[order[2]])
        lb = np.sqrt(D[order[1]])
        lc = np.sqrt(D[order[0]])

        dq = np.abs(q-lb/la)
        ds = np.abs(s-lc/la)

        q = lb/la
        s = lc/la

        Tiv = np.dot(LA.inv(A), Tiv)

    # rba = q
    # rca = s
    if decrease:
        Tiv = Tiv[order[::-1], :]
    else:
        Tiv = Tiv[order, :]
    # right hand coordinate
    if np.dot(Tiv[2, :], np.cross(Tiv[0, :], Tiv[1, :])) < 0:
        Tiv[1, :] *= -1
    # eigen vectors Vei[:,0] Vei[:,1] Vei[:,2]
    Vei = LA.inv(Tiv)

    d = np.array([0, 0, 1])
    costh = np.dot(Vei[:, 2], d) / np.sqrt(np.dot(Vei[:, 2], Vei[:, 2])) /\
        np.sqrt(np.dot(d, d))
    # angle between longest axis (z' direction) and LOS (i.e. z direction)
    angle = np.arccos(costh)*180./np.pi

    ba = q
    ca = s
    # print ' b/a= {:.2f}'.format(q),'  c/a = {:.2f}'.format(s)
    # print "rotation angle=",angle

    return ba, ca, angle, Tiv


def http_get(path, params=None):
    # make HTTP GET request to path
    headers = {"api-key": "c7c3eb002b50d82782fcfc277043074e"}
    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json()  # parse json responses automatically

    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename  # return the filename string

    return r


def http_info2dict(info):
    dictionary = {}
    for key in info.keys():
        dictionary[key] = info[key]
    return dictionary
