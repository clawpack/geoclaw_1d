

import os, sys
from imp import reload


try:
    from clawpack.geoclaw_1d import geoplot
except:
    print('Could not import from geoclaw_1d')

#import clawpack.geoclaw.shallow_1d.plot as geoplot

import setrun
rundata=setrun.setrun()

import mapc2p
reload(mapc2p)  # in case num_cells changed
from mapc2p import make_mapc2p

import numpy
#from pylab import find

try:
    fname = '_output/fort.hmax'
    d = numpy.loadtxt(fname)
    etamax = numpy.where(d[:,1]>1e-6, d[:,2], numpy.nan)
    xmax = d[:,0]
    jmax = where(d[:,1]>0)[0].max()
    print("run-in = %8.2f,  run-up = %8.2f" % (d[jmax,0],d[jmax,2]))
    print('Loaded hmax from ',fname)
except:
    xmax = None
    print("Failed to load fort.hmax")

xlimits = [-150e3,150e3]

def setplot(plotdata):

    plotdata.clearfigures()

    outdir1 = plotdata.outdir
    mapc2p1, ngrid1 = make_mapc2p(outdir1)


    def fixticks1(current_data):
        from pylab import ticklabel_format, grid,tight_layout
        ticklabel_format(useOffset=False)
        grid(True)
        tight_layout()

    def fixticks(current_data):
        from pylab import ticklabel_format, plot,grid,gca
        ticklabel_format(useOffset=False)
        if xmax is not None:
            plot(xmax, etamax, 'r')
        grid(True)

    plotfigure = plotdata.new_plotfigure(name='domain', figno=0)
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = xlimits
    #plotaxes.xlimits = [-100e3,-20e3]
    plotaxes.ylimits = [-1,4]
    plotaxes.title = 'Surface displacement'
    plotaxes.afteraxes = fixticks
    plotaxes.skip_patches_outside_xylimits = False

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.surface
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'k'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = xlimits
    #plotaxes.xlimits = [-100e3,-20e3]
    #plotaxes.ylimits = [-1000, 1000]
    #plotaxes.title = 'Full depth'
    plotaxes.title = 'momentum'
    plotaxes.afteraxes = fixticks1
    plotaxes.skip_patches_outside_xylimits = False
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.show = False
    plotitem.plot_var = geoplot.surface
    plotitem.plot_var2 = geoplot.topo
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'k'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 1
    plotitem.color = 'k'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    #----------

    plotfigure = plotdata.new_plotfigure(name='shore', figno=1)
    #plotfigure.kwargs = {'figsize':(9,11)}
    plotfigure.show = False
    

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = [0,80e3]
    plotaxes.ylimits = [-4,4]
    plotaxes.title = 'Zoom on shelf'

    plotaxes.afteraxes = fixticks
    plotaxes.skip_patches_outside_xylimits = False

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.surface
    #plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    #plotitem.plot_var = geoplot.surface
    #plotitem.plot_var2 = geoplot.topo
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'k'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    #plotaxes.xlimits = [-2000,2000]
    plotaxes.xlimits = [-1000,1000]
    #plotaxes.ylimits = [-10,40]
    plotaxes.ylimits = [-20,60]
    plotaxes.title = 'Zoom around shore'

    plotaxes.afteraxes = fixticks
    plotaxes.skip_patches_outside_xylimits = False

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False
    plotitem.plot_var = geoplot.surface

    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.plot_var = geoplot.surface
    plotitem.plot_var2 = geoplot.topo
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'k'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p



    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='q', figno=300, \
                                         type='each_gauge')
    plotfigure.clf_each_gauge = True

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-2,2]
    plotaxes.title = 'Surface elevation eta'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 2  # eta
    plotitem.plotstyle = 'b-'


    plotdata.printfigs = True          # Whether to output figures
    plotdata.print_format = 'png'      # What type of output format
    plotdata.print_framenos = 'all'      # Which frames to output
    plotdata.print_fignos = 'all'      # Which figures to print
    plotdata.html = True               # Whether to create HTML files
    plotdata.latex = False             # Whether to make LaTeX output
    plotdata.parallel = True

    return plotdata

