
.. _geoclaw_1d/examples/okada_dtopo:

Earthquake generated tsunami using the Okada model
===================================================================


The file `make_celledges.py` sets up the domain, computational grid, and
topography, creating a file celledges.data.

A piecewise linear topography is defined by specifying the topography `z`
value at a set of nodes `x` in the `xzpairs` list.  
In this example is is set up to give a 4000 m deep ocean with a continental
slope leading to a 200 m deep shelf, followed by a beach, as shown in
topo.png after running make_celledges.py.

A nonuniform grid with `mx` grid cells is used with cell widths related
to the still water depth in such a way that the Courant number is roughly
constant in deep water and onto the shelf, and with uniform grid cells
near shore and onshore where the water depth is less than `hmin`.

Executing `make_celledges.py` creates a file `celledges.data` that contains
the cell edges and also topography values at these points.
This file must be created before running GeoClaw.

Note that this file is used also as the topofile, as specified in setrun.py.

In GeoClaw a mapped grid is used with a `mapc2p` function specified in
`setrun.py` that is generated from the `celledges.data`.  The computational
grid specified in `setrun.py` is always `0 <= xc <= 1`.  Set::

    rundata.grid_data.grid_type = 2
    
to indicate a mapped grid.

In this example the physical `x` coordiate is in meters, set by specifying::

    rundata.geo_data.coordinate_system = 1

This example also uses a dtopo file generated by the code make_dtopo.py.
The Okada model is used to compute the dtopo on a subducting fault, with the
one-dimensional direction modeled here corresponding to the direction of the
fault dip.  This models a fault that is infinitely long in the strike
direction orthogonal to this 1d model, also assuming the topography is
constant in that direction.

Executing make_dtopo.py also produces plots of the subfaults and seafloor
deformation.

To use::

    make topo     # executes make_celledges.py
    make dtopo     # executes make_dtopo.py
    make .output  # compile, make data, and run
    make .plots   # to create _plots (or plot interactively with Iplotclaw)


Version
-------

July 2023, For initial release of geoclaw 1d code in Clawpack v5.10
