
from __future__ import print_function
import numpy as np
import os
from scipy.interpolate import interp1d

def make_mapc2p(outdir):
    """
    Create a mapc2p function that maps computational cell edges xc
    with 0 <= xc <= 1 to the physical cell edges.  The physical
    cell edges should be in the file grid.data, starting in the third row.
    """

    grid_data_file = os.path.join(outdir, 'grid.data')
    d = np.loadtxt(grid_data_file, skiprows=3) 
    ngrid = d.shape[0]

    print('mapc2p: Read %i grid values from %s' % (d.shape[0], grid_data_file))

    xc_edges = np.linspace(0,1,ngrid)
    xp_edges = d[:,0]
    mapc2p_fcn = interp1d(xc_edges, xp_edges, kind='linear')

    return mapc2p_fcn, ngrid
