"""
Plots time-dependent Okada solution
"""

from pylab import *

import setfault
reload(setfault)

fault = setfault.make_fault()

tend = 0.
for s in fault.subfaults:
    tend = max(tend, s.rupture_time + 2*s.rise_time)
times = linspace(0, tend+1., 10)

x = linspace(-2,2,1000)
y = array([0,1])

dtopo = fault.create_dtopography(x,y,times)


figure(350)
clf()
for k in range(len(times)):
    plot(111e3*dtopo.X[0,:],dtopo.dZ[k,0,:],label='t = %6.2f' % times[k])
    
legend(loc='upper right')
title('Okada solution')
