"""
Plot results at three gauges and compare with experimental data, obtained from
    https://nctr.pmel.noaa.gov/benchmark/Solitary_wave/
Compare this figure to Figure 4 in the BoussClaw paper
    http://dx.doi.org/10.1016/j.coastaleng.2017.01.005
"""

from pylab import *
import clawpack.pyclaw.gauges as gauges

d = loadtxt('experimental_data/ts3b.txt',skiprows=6)
tg = d[:,0]
g5 = d[:,2]
g7 = d[:,4]
g8 = d[:,5]


figure(400, figsize=(8,8))
clf()

subplot(311)
gaugeno = 5
gauge = gauges.GaugeSolution(gaugeno, '_output')
t = gauge.t
eta = gauge.q[2,:]

plot(tg-tg[0], g5, 'r', label='Experiment')
plot(t, eta, 'b', label='Bouss')
xlim(0,25)
ylim(-0.01,0.08)
legend(loc='upper right')
grid(True)
xlabel('')
ylabel('Surface (m)')
title('Gauge %i' % gaugeno)

subplot(312)
gaugeno = 7
gauge = gauges.GaugeSolution(gaugeno, '_output')
t = gauge.t
eta = gauge.q[2,:]

plot(tg-tg[0], g7, 'r', label='Experiment')
plot(t, eta, 'b', label='Bouss')
xlim(0,25)
ylim(-0.01,0.08)
legend(loc='upper right')
grid(True)
xlabel('')
ylabel('Surface (m)')
title('Gauge %i' % gaugeno)

subplot(313)
gaugeno = 8
gauge = gauges.GaugeSolution(gaugeno, '_output')
t = gauge.t
eta = gauge.q[2,:]

plot(tg-tg[0], g8, 'r', label='Experiment')
plot(t, eta, 'b', label='Bouss')
xlim(0,25)
ylim(-0.01,0.08)
legend(loc='upper right')
grid(True)
xlabel('')
ylabel('Surface (m)')
title('Gauge %i' % gaugeno)

tight_layout()


if 1:
    fname = 'GaugeComparison.png'
    savefig(fname, bbox_inches='tight')
    print('Created %s' % fname)

