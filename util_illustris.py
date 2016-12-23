import numpy as np
import h5py
import numpy.linalg as LA

# global parameters
boxsize_img = 50.0
scale_img = 0.25  # pixel2kpc
boxsize_ifu = 50.0
scale_ifu = 0.31
kpc2arcsec = 1.612


def read_cutout(fname, PartType=4, key='Coordinates', z=0.0):
    '''
    keys:
      dm - Coordinates, Velocities
      stars - Coordinates, Velocities, GFM_Metallicity,
              GFM_StellarFormationTime, GFM_StellarPhotometrics, Masses,
    '''
    little_h = 0.704
    scale_factor = 1.0 / (1+z)
    with h5py.File(fname, 'r') as f:
        Ptype = 'PartType{}'.format(PartType)
        data = f[Ptype][key][:]
    if key == 'Coordinates':
        # distance in kpc
        data *= scale_factor/little_h
    elif key == 'Velocities':
        # velocity in km/s
        data *= scale_factor
    elif key == 'Masses':
        # mass in M_sun
        data /= little_h
    return data


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
