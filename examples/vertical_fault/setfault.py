from __future__ import print_function
from numpy import arange,cos,sin,pi
from clawpack.geoclaw.data import LAT2METER

import clawpack.seismic.dtopotools_horiz_okada_and_1d as dtopotools
#reload(dtopotools)

def make_fault():
    fault = dtopotools.Fault(coordinate_specification='top center')
    fault.subfaults = []

    width = 10e3
    #theta = 0.20
    dip = 90.  # dip in degrees
    fault_top_center = [0.,-10e3]
    average_slip = 2.0
    max_slip = 2*average_slip # if modulated by cosine hump below
    mu = 3e10
    rupture_time = 0.0
    rise_time = 20.
    nsubfaults = 20

    theta = dip*pi/180.
    longitude0 = fault_top_center[0]/LAT2METER
    dlongitude = width*cos(theta)/LAT2METER / nsubfaults
    ddepth = width*sin(theta) / nsubfaults
    subfault_width = width/nsubfaults

    total_slip = 0.
    for i in range(nsubfaults):
        # split total slip between subfaults, starting at top
        subfault = dtopotools.SubFault()
        subfault.mu = mu
        subfault.dip = dip  #theta/pi*180.0
        subfault.width = subfault_width
        subfault.depth = -fault_top_center[1] + ddepth*i
        subfault.slip = max_slip * 0.5*(1 - cos(2*pi*(i+0.5)/nsubfaults))
        #subfault.slip = average_slip
        total_slip += subfault.slip
        print('subfault %2i at depth %8.3f km has slip = %6.3f' \
                % (i,subfault.depth/1e3,subfault.slip))
        subfault.rake = 90
        subfault.strike = 0
        subfault.length = 1000e3
        subfault.longitude = longitude0 + dlongitude*i
        subfault.latitude = 0.
        subfault.coordinate_specification = 'top center'
        subfault.rupture_time = rupture_time
        subfault.rise_time = rise_time
        fault.subfaults.append(subfault)

    print('average slip = %6.3f' % (total_slip/nsubfaults))
    fault.rupture_type = 'dynamic'
    return fault

if __name__=='__main__':
    fault = make_fault()
    fault.write('fault.data')
