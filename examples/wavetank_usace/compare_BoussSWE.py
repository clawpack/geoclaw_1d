"""
Run GeoClaw with both Bouss terms and pure shallow water equations,
and then plot gauge results compared to experiment.

"""
from pylab import *
import os

from clawpack.clawutil.runclaw import runclaw
#from clawpack.visclaw.frametools import plotframe
#from clawpack.visclaw.data import ClawPlotData
import compare_gauges

import setrun

outdir_bouss = '_output_bouss'
print('outdir_bouss = ',outdir_bouss)

outdir_swe = '_output_swe'
print('outdir_swe = ',outdir_swe)

run_code = False # set to False if output already exists

if run_code:
    # create executable and .data files:
    os.system('make .exe')
    rundata = setrun.setrun()

    # Boussinesq:
    rundata.bouss_data.bouss = True
    rundata.write()
    runclaw(xclawcmd='xgeo',outdir=outdir_bouss)   # run clawpack code
    
    # Shallow water equations:
    rundata.bouss_data.bouss = False
    rundata.write()
    runclaw(xclawcmd='xgeo',outdir=outdir_swe)   # run clawpack code
    

savefig_ext = '.png'
figdir = './figures'
os.system('mkdir -p %s' % figdir)

def save_figure(fname):
    """Save figure to figdir with desired extension"""
    full_fname = os.path.join(figdir,fname) + savefig_ext
    savefig(full_fname, bbox_inches='tight')
    print('Created %s' % full_fname)

compare_gauges.plot_gauges(outdir_bouss, outdir_swe,
                           'GaugeComparison_BoussSWE.png')
