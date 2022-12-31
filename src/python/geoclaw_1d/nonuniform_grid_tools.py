from pylab import *
from scipy.interpolate import interp1d
import numpy as np
import os


def make_mapc2p(fname_celledges='celledges.data'):
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

    xc_edge = np.linspace(0,1,mx_edge)  # assumes xlower=0, xupper=1 in setrun
    xp_edge = d[:,0]
    mapc2p = interp1d(xc_edge, xp_edge, kind='linear')

    return mapc2p, mx_edge, xp_edge


def make_pwlin_topo_fcn(xzpairs):
    xi = array([xz[0] for xz in xzpairs])
    zi = array([xz[1] for xz in xzpairs])
    z_fcn = interp1d(xi, zi, kind='linear', bounds_error=False, 
                     fill_value='extrapolate')
    return z_fcn


def make_celledges_cfl(xlower, xupper, mx, topo_fcn, hmin,
                       fname='celledges.data', plot_topo=False):

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
    dxp = diff(xp)
    
    if plot_topo:
        figure(97, figsize=(8,8))
        clf()
        subplot(211)
        #plot(csum, xunif, 'b')
        plot(xunif, csum, 'b')
        ylabel('computational coordinate xc')
        grid(True)
        axis([xlower,xupper,0,1])
        title('inverse of mapc2p function')

        subplot(212)
        xcell = 0.5*(xp[1:] + xp[:-1])
        plot(xcell, dxp, 'b')
        xlabel('physical coordinate xp')
        ylabel('delta x')
        grid(True)
        xlim(xlower,xupper)
        title('Mesh width')
        tight_layout()
        
        png_fname = 'cellmap.png'
        savefig(png_fname)
        print("Created ",png_fname)

    with open(fname,'w') as f:
        f.write('%i   # number of cell edges\n' % (mx+1))

        for i in range(mx+1):
            f.write('%15.8f %15.8f\n' % (xp[i],z[i]))
        f.close()

    print("Created %s, containing %i cell edges" % (fname,mx+1))
    print("Min dx = %g, Max dx = %g" % (dxp.min(),dxp.max()))

    if plot_topo:
        figure(99, figsize=(8,4))
        clf()
        fill_between(xp,where(z<0,z,nan),0.,color=[.5,.5,1])
        plot(xp,z,'g')
        xlim(xlower,xupper)
        zmax = max(z.max(), 0)
        zmargin = 0.1*(zmax-z.min())
        ylim(z.min()-zmargin,zmax+zmargin)
        grid(True)
        title('Topography')
        png_fname = 'topo.png'
        savefig(png_fname)
        print("Created ",png_fname)

