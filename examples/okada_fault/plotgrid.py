
from pylab import *

figure(1); clf()
griddata = loadtxt('grid.data', skiprows=1)
xgrid = griddata[:,0]
zgrid = griddata[:,1]
plot(xgrid,zgrid,'g')
ylim(-4500,0)
xlim(-150e3,150e3)

