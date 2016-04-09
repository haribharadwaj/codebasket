# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 14:36:30 2015

@author: Hari
"""
import pyNN.neuron as sim
import pylab as pl
import numpy as np
from pyNN.random import RandomDistribution


sim.setup(timestep=0.01)
nE = 40
nI = 20
Ecells = sim.Population(nE, sim.HH_cond_exp(tau_syn_E=1.0, tau_syn_I=5.0),
                        label='E-Cells')
Icells = sim.Population(nI, sim.HH_cond_exp(tau_syn_E=1.0, tau_syn_I=5.0),
                        label='I-Cells')
Net = Ecells + Icells
vinit_distr = RandomDistribution('uniform', (-70., -55.))
Net.initialize(v=vinit_distr)  # Initialize to random subthreshold voltages

Esyn = sim.StaticSynapse(weight=0.05, delay=1.0)
random = sim.FixedProbabilityConnector(p_connect=0.3)
Econnections = sim.Projection(Ecells, Net, random, Esyn,
                              receptor_type='excitatory')
Isyn = sim.StaticSynapse(weight=0.10, delay=1.0)
Iconnections = sim.Projection(Icells, Net, random, Isyn,
                              receptor_type='inhibitory')

noise1 = sim.NoisyCurrentSource(mean=0.6, stdev=0.3, start=50.0,
                                stop=850.0, dt=0.01)

Ecells.inject(noise1)
# noise2 = sim.NoisyCurrentSource(mean=0.3, stdev=0.3, start=50.0,
#                                 stop=850.0, dt=0.01)

# noise1.inject_into(Ecells[(0, 5, 10, 11)])
# noise2.inject_into(Ecells[(3, 7, 15, 17)])

# sine = sim.ACSource(start=50.0, stop=850.0, amplitude=1.0, offset=1.0,
#                     frequency=25.0, phase=180.0)

# sine.inject_into(Ecells[(2, 25, 39)])

Net.record(['v', 'spikes'])
sim.run(1000.0)

data = Net.get_data()  # Neo data format for electrophysiology data

sim.end()

# Plot Membrane Voltages
V = data.segments[0].analogsignalarrays[0]
# Plot E-cells
ax1 = pl.subplot(211)
for k in range(nE):
    pl.plot(V.times, V[:, k], linewidth=2)
    pl.hold(True)

pl.xlabel('Time (ms)', fontsize=16)
pl.ylabel('Membrane potential (mV)', fontsize=16)
pl.title('E-Cells', fontsize=16)
pl.tick_params(labelsize=16)


# Plot I-cells
ax2 = pl.subplot(212, sharex=ax1)
for k in range(nE, nE + nI):
    pl.plot(V.times, V[:, k], linewidth=2)
    pl.hold(True)

pl.xlabel('Time (ms)', fontsize=16)
pl.ylabel('Membrane potential (mV)', fontsize=16)
pl.title('I-Cells', fontsize=16)
pl.tick_params(labelsize=16)
pl.show()

# Plot Spike Times
spikes = data.segments[0].spiketrains
pl.figure()
for k in range(nE + nI):
    spike = spikes[k]
    col = '.k' if k < nE else '.r'
    if k == 0:
        lab = 'E-Cells'
    else:
        if k == nE:
            lab = 'I-Cells'
        else:
            lab = None
    pl.plot(spike, np.ones_like(spike) * k, col, markersize=200/nE,
            label=lab)
    pl.hold(True)

pl.xlabel('Time (ms)', fontsize=16)
pl.ylabel('Neuron Index', fontsize=16)
pl.title('Simple E-I network', fontsize=16)
pl.tick_params(labelsize=16)
pl.legend(loc='best')
pl.show()
