"""
Run GeoClaw with pure shallow water equations and with several choices
of Boussinesq equations, and then plot gauge results compared to experiment.

"""
from pylab import *
import os

from clawpack.clawutil.runclaw import runclaw
#from clawpack.visclaw.frametools import plotframe
#from clawpack.visclaw.data import ClawPlotData
import compare_gauges

import setrun

outdir_sgn = '_output_sgn'
print('outdir_bouss = ',outdir_sgn)

outdir_ms = '_output_ms'
print('outdir_ms = ',outdir_ms)

outdir_swe = '_output_swe'
print('outdir_swe = ',outdir_swe)

run_code = True  # set to False if output already exists

if run_code:
    # create executable and .data files:
    os.system('make .exe')
    rundata = setrun.setrun()

    # Boussinesq, MS:
    rundata.bouss_data.boussEquations = 1
    rundata.write()
    runclaw(xclawcmd='xgeo',outdir=outdir_ms)   # run clawpack code

    # Boussinesq, SGN:
    rundata.bouss_data.boussEquations = 2
    rundata.write()
    runclaw(xclawcmd='xgeo',outdir=outdir_sgn)   # run clawpack code
    
    # Shallow water equations:
    rundata.bouss_data.boussEquations = 0
    rundata.write()
    runclaw(xclawcmd='xgeo',outdir=outdir_swe)   # run clawpack code
    

savefig_ext = '.pdf'
figdir = './figures'
os.system('mkdir -p %s' % figdir)

def save_figure(fname):
    """Save figure to figdir with desired extension"""
    full_fname = os.path.join(figdir,fname) + savefig_ext
    savefig(full_fname, bbox_inches='tight')
    print('Created %s' % full_fname)

outdirs=[('_output_swe', 'SWE', 'k'), \
         ('_output_ms', 'MS, B = 1/15', 'b'), \
         ('_output_sgn', 'SGN, alpha = 1.153','g')]

compare_gauges.plot_gauges(outdirs, fname_figure='GaugeComparison_BoussSWE.png')
