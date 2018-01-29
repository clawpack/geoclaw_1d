
# okada fault

Initial test using Okada model to generate initial displacement
used in qinit (we don't yet have moving topography incorporated).

Uses a nonuniform grid designed so that the Courant number is close to 1
everywhere.  Create grid.data by:

    python makegrid.py    # or:  make grid

This also calls make_qinit_okada.py to create qinit_okada.data, the Okada
solution evaluated on this grid.  The Okada parameters are set in 
setfault.py.


