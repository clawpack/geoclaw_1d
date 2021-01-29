from pylab import *
from scipy.interpolate import interp1d
from interp import pwcubic, pwlinear
import numpy as np
import os
from scipy.interpolate import interp1d


def make_mapc2p(fname_celledges='celledges.txt'):
    """
    Create a mapc2p function that maps computational cell edges xc
    with 0 <= xc <= 1 to the physical cell edges.  The physical
    cell edges should be in the file fname_celledges, 
    starting in the second row (following mx_edge, the number of edges).

    Returns the mapc2p function and also mx_edge, xp_edge, which may be
    needed for other purposes.
    """

    path = os.path.abspath(fname_celledges)
    d = np.loadtxt(path, skiprows=1) 
    mx_edge = d.shape[0]

    print('make_mapc2p: Using %i cell edge values from %s' % (d.shape[0], path))

    xc_edge = np.linspace(0,1,mx_edge)
    xp_edge = d[:,1]
    mapc2p = interp1d(xc_edge, xp_edge, kind='linear')

    return mapc2p, mx_edge, xp_edge


def make_pwlin_topo_fcn(xzpairs):
    xi = array([xz[0] for xz in xzpairs])
    zi = array([xz[1] for xz in xzpairs])
    z_fcn = interp1d(xi, zi, kind='linear', bounds_error=False, 
                     fill_value='extrapolate')
    return z_fcn


def make_celledges_cfl(xlower, xupper, mx, topo_fcn, hmin,
                       cfl=1., fname='celledges.txt', plot_topo=False):

    grav = 9.81

    cmin = sqrt(grav*hmin)

    def c(x):
        z = topo_fcn(x)
        h = where(-z > hmin, -z, hmin)
        c = sqrt(grav*h)
        return c

    xunif = linspace(xlower, xupper, 2*mx)
    cunif = c(xunif)
    csum = cumsum(1./cunif)
    csum = csum - csum[0]

    csum = csum / csum[-1]
    cinv = interp1d(csum, xunif)

    xc = linspace(0, 1, mx+1)   # computational grid
    xp = cinv(xc)
    z = topo_fcn(xp)

    with open(fname,'w') as f:
        f.write('%i   # number of cell edges\n' % (mx+1))

        for i in range(mx+1):
            f.write('%15.4f %15.4f %15.4f\n' % (xc[i],xp[i],z[i]))
        f.close()

    print("Created %s, containing cell edges" % fname)

    if plot_topo:
        figure(1, figsize=(8,4))
        clf()
        fill_between(xp,where(z<0,z,nan),0.,color=[.5,.5,1])
        plot(xp,z,'g')
        xlim(xlower,xupper)
        ylim(z.min()-500,500)
        fname = 'topo.png'
        savefig(fname)
        print("Created ",fname)

