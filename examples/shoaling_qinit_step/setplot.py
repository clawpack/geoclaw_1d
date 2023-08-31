

import os, sys

try:
    from clawpack.geoclaw_1d import geoplot
except:
    print('Could not import from geoclaw_1d')


from clawpack.geoclaw_1d.nonuniform_grid_tools import make_mapc2p
import numpy

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

fname_celledges = os.path.abspath('celledges.data')


xlimits = [-300,50]

def setplot(plotdata=None):

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()

    plotdata.clearfigures()

    mapc2p, mx_edge, xp_edge = make_mapc2p(fname_celledges)

    def mapc2p_km(xc):
        x_m = mapc2p(xc)
        x_km = x_m / 1000.   # convert to km
        return x_km


    def fixticks1(current_data):
        from pylab import ticklabel_format, grid
        ticklabel_format(useOffset=False)
        grid(True)

    def fixticks(current_data):
        from pylab import ticklabel_format, plot,grid,ones,sqrt, \
            legend,title,ylabel,text
        ticklabel_format(useOffset=False)

        # to plot max elevation over entire computation:
        #if xmax is not None:
        #    plot(xmax, etamax, 'r')

        #grid(True)
        hl = 3200.
        hr = 200.
        greens = (hl/hr)**(0.25)
        print('greens = ',greens)
        #plot(current_data.x, greens*ones(current_data.x.shape),'g--')
        plot(xlimits,[greens,greens],'g--', label='$C_g$, Greens Law')
        ctrans = 2*sqrt(hl)/(sqrt(hl)+sqrt(hr))
        crefl = (sqrt(hl)-sqrt(hr))/(sqrt(hl)+sqrt(hr))
        print('ctrans = ',ctrans)
        plot(xlimits,[ctrans,ctrans],'r--', label='$C_T$, Transmission coefficient')
        legend(loc='upper left')
        title('')
        ylabel('meters', fontsize=14)

        if current_data.frameno == 0:
            text(-80,-0.4,'$\longrightarrow$',fontsize=20)
            text(-80,-0.6,'Incident')
        h = current_data.q[0,:]
        mx2 = int(round(len(h)/2.))
        etamax2 = (h[:mx2] - hl).max()
        print('mx2 = %i, etamax2 = %g' % (mx2,etamax2))

    plotfigure = plotdata.new_plotfigure(name='domain', figno=0)
    plotfigure.kwargs = {'figsize':(7,6.5)}
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'axes([.1,.4,.8,.5])' #'subplot(211)'
    plotaxes.xlimits = xlimits
    #plotaxes.xlimits = [-100e3,-20e3]
    plotaxes.ylimits = [-1,3]
    plotaxes.title = 'Surface displacement'
    plotaxes.afteraxes = fixticks

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.surface
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p_km


    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False
    plotitem.plot_var = 1
    plotitem.color = 'k'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p_km

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'axes([.1,.1,.8,.2])' #'subplot(212)'
    plotaxes.xlimits = xlimits

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'g'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p_km

    def fix_topo_plot(current_data):
        from pylab import title,xlabel
        title('')
        xlabel('kilometers', fontsize=14)
    plotaxes.afteraxes = fix_topo_plot

    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p_km



    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='q', figno=300, \
                                         type='each_gauge')
    plotfigure.clf_each_gauge = True

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Eta'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 2
    plotitem.plotstyle = 'b-'

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:


    plotdata.printfigs = True          # Whether to output figures
    plotdata.print_format = 'png'      # What type of output format
    plotdata.print_framenos = 'all'      # Which frames to output
    plotdata.print_fignos = 'all'      # Which figures to print
    plotdata.html = True               # Whether to create HTML files
    plotdata.latex = False             # Whether to make LaTeX output
    plotdata.parallel = True

    return plotdata

