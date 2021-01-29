
from pylab import *
from clawpack.geoclaw_1d import nonuniform_grid_tools


xzpairs = [(-150e3,-4500),   # left edge
           (  50e3,-4500),   # start of continental slope
           ( 100e3,-480),    # start of continental shelf
           ( 150e3,-480)]    # right edge

topo_fcn = nonuniform_grid_tools.make_pwlin_topo_fcn(xzpairs)

mx = 5000

nonuniform_grid_tools.make_celledges_cfl(-150e3,150e3,mx,topo_fcn,
        hmin=50, cfl=1, fname='celledges.txt', plot_topo=True)

