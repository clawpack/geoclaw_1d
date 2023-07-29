
from numpy import arange,cos,sin,pi
from clawpack.geoclaw.data import LAT2METER
from importlib import reload

from clawpack.geoclaw_1d import dtopotools
reload(dtopotools)

def make_fault():
    fault = dtopotools.Fault(coordinate_specification='top center')
    fault.subfaults = []

    fault_top_meters = -100e3
    fault_top_depth = 20e3
    width = 50e3

    theta = 0.20  # dip in radians
    dip = theta/pi*180.0

    #dip = 15.  # dip in degrees
    #theta = dip*pi/180.

    average_slip = 1.0
    max_slip = 2*average_slip # if modulated by cosine hump below
    mu = 3e10
    rupture_time = 0.0
    rise_time = 10.
    nsubfaults = 2

    longitude0 = fault_top_meters/LAT2METER
    dlongitude = width*cos(theta)/LAT2METER / nsubfaults
    ddepth = width*sin(theta) / nsubfaults
    subfault_width = width/nsubfaults

    total_slip = 0.
    for i in range(nsubfaults):
        # split total slip between subfaults, starting at top
        subfault = dtopotools.SubFault()
        subfault.mu = mu
        subfault.dip = dip
        subfault.width = subfault_width
        subfault.depth = fault_top_depth + ddepth*i
        #subfault.slip = max_slip * 0.5*(1 - cos(2*pi*(i+0.5)/nsubfaults))
        subfault.slip = average_slip  # for constant slip
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
