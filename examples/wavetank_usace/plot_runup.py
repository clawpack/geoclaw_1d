
from pylab import *

d = loadtxt('_output/runup.txt')
t = d[:,0]
xshore = d[:,3]
zshore = d[:,4]

figure(80,figsize=(10,6))
clf()
plot(t, zshore, 'b')
plot(t, xshore, 'g')
grid(True)
title('Runup z vs. time')
