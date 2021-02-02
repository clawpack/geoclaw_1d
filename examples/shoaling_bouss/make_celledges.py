
from pylab import *
from clawpack.geoclaw_1d import nonuniform_grid_tools


xzpairs = [(   0e3,-1000),   # left edge
           (  20e3,-1000),   # start of continental slope
           (  25e3,-200),    # start of continental shelf
           (  50e3,-200),    # start of beach
           (  60e3,0),       # shore
           (  64e3,80)]      # right edge

topo_fcn = nonuniform_grid_tools.make_pwlin_topo_fcn(xzpairs)

mx = 1000

nonuniform_grid_tools.make_celledges_cfl(0, 64e3, mx, topo_fcn,
        hmin=50, cfl=1, fname='celledges.txt', plot_topo=True)

