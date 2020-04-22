"""
Module to set up run time parameters for Clawpack -- classic code.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""

from __future__ import print_function

import os, sys
import numpy as np
from mapc2p import mapc2p

dispersion = True            # Include Boussinesq terms?
B_param =  1.0 / 15.0         # Parameter for the Boussinesq eqns
sw_depth0 = 20.              # Use pure SWE if depth less than sw_depth0
sw_depth1 = 20.              # Use pure Bous if depth greater than sw_depth1
radial = True               # Include radial source terms?

grid = np.loadtxt('grid.data', skiprows=1)
print('Read grid from grid.data, %i grid values' % grid.shape[0])
mx = grid.shape[0] - 1

dxc = 1./mx
xc = np.linspace(dxc/2., 1-dxc/2., mx)
xp = mapc2p(xc)

#------------------------------
def setrun(claw_pkg='geoclaw'):
#------------------------------

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    """

    from clawpack.clawutil import data

    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    num_dim = 1
    rundata = data.ClawRunData(claw_pkg, num_dim)

    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------
    # Sample setup to write one line to setprob.data ...
    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')
    probdata.add_param('bouss'             , dispersion  ,'Include dispersive terms?')
    probdata.add_param('B_param'           , B_param     ,'Parameter for the Boussinesq eq')
    probdata.add_param('sw_depth0', sw_depth0)
    probdata.add_param('sw_depth1', sw_depth1)
    probdata.add_param('radial', radial, 'Radial source term?')



    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #------------------------------------------------------------------

    clawdata = rundata.clawdata  # initialized when rundata instantiated


    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim

    # Lower and upper edge of computational domain:
    clawdata.lower[0] = 0.          # xlower
    clawdata.upper[0] = 1.           # xupper

    # Number of grid cells:
    clawdata.num_cells[0] = mx


    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 2

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 2

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 2


    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.


    # Restart from checkpoint file of a previous run?
    # Note: If restarting, you must also change the Makefile to set:
    #    RESTART = True
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.qNNNN' specified below should be in
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False               # True to restart from prior results
    clawdata.restart_file = 'fort.q0006'   # File to use for restart data


    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.

    clawdata.output_style = 1

    if clawdata.output_style==1:
        # Output ntimes frames at equally spaced times up to tfinal:
        # Can specify num_output_times = 0 for no output
        clawdata.num_output_times = 30
        clawdata.tfinal = 3000.
        clawdata.output_t0 = False  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list or numpy array of output times:
        # Include t0 if you want output at the initial time.
        clawdata.output_times =  [0.2,1.6,2.9,4.3,5.6,6.9,8.3,9.6,10.9]

    elif clawdata.output_style == 3:
        # Output every step_interval timesteps over total_steps timesteps:
        clawdata.output_step_interval = 10
        clawdata.total_steps = 100
        clawdata.output_t0 = True  # output at initial (or restart) time?


    clawdata.output_format = 'ascii'      # 'ascii', 'binary', 'netcdf'

    clawdata.output_q_components = 'all'   # could be list such as [True,True]
    clawdata.output_aux_components = 'all'  # could be list
    clawdata.output_aux_onlyonce = True   # output aux arrays only at t0


    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 0



    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==True:  variable time steps used based on cfl_desired,
    # if dt_variable==False: fixed time steps dt = dt_initial always used.
    clawdata.dt_variable = True

    # Initial time step for variable dt.
    # (If dt_variable==0 then dt=dt_initial for all steps)
    clawdata.dt_initial = 0.01

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1.e9

    # Desired Courant number if variable dt used
    clawdata.cfl_desired = 0.75
    # max Courant number to allow without retaking step with a smaller dt:
    clawdata.cfl_max = 1.0

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 50000


    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2


    # Number of waves in the Riemann solution:
    clawdata.num_waves = 2

    # List of limiters to use for each wave family:
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'vanleer'  ==> van Leer
    #   4 or 'mc'       ==> MC limiter
    clawdata.limiter = [4,4]

    clawdata.use_fwaves = True    # True ==> use f-wave version of algorithms

    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used,
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 'godunov'


    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    # Choice of BCs at xlower and xupper:
    #   0 or 'user'     => user specified (must modify bcNamr.f to use this option)
    #   1 or 'extrap'   => extrapolation (non-reflecting outflow)
    #   2 or 'periodic' => periodic (must specify this at both boundaries)
    #   3 or 'wall'     => solid wall for systems where q(2) is normal velocity

    #clawdata.bc_lower[0] = 'extrap'   # at xlower
    clawdata.bc_lower[0] = 'wall'   # at xlower
    clawdata.bc_upper[0] = 'extrap'   # at xupper


    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft'  (see documentation).
    # Isn't used for this non-amr version, but still expected in data.

    amrdata = rundata.amrdata
    amrdata.aux_type = ['center','capacity']


    # ---------------
    # Gauges:
    # ---------------
    rundata.gaugedata.gauges = []
    # for gauges append lines of the form  [gaugeno, x, t1, t2]
    xgauge = [-50e3, 0, 50e3, 100e3]
    for gaugeno,xp_g in enumerate(xgauge):
        # compute computational point xc_g that maps to xp_g:
        ii = np.where(xp < xp_g)[0][-1]
        xp_frac = (xp_g - xp[ii])/(xp[ii+1] - xp[ii])
        xc_g = (ii + xp_frac)/float(grid.shape[0])
        print('gaugeno = ',gaugeno)
        print('  ii, xp_g, xp[ii], xp[ii+1], xp_frac, xc_g:\n ', \
                 ii, xp_g, xp[ii], xp[ii+1], xp_frac, xc_g)
        rundata.gaugedata.gauges.append([gaugeno, xc_g, 0, 1e9])


    # ---------------
    # geo data
    # ---------------

    geo_data = rundata.geo_data

    geo_data.sea_level = 0.0
    geo_data.dry_tolerance = 1.e-3

    # Friction source terms:
    #   src_split > 0 required

    geo_data.friction_forcing = True
    geo_data.manning_coefficient =.025
    geo_data.friction_depth = 1e6


    return rundata

    # end of function setrun
    # ----------------------


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    rundata = setrun(*sys.argv[1:])
    rundata.write()

