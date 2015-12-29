# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 14:36:30 2015

@author: Hari
"""
import pyNN.neuron as sim
import pylab as pl

sim.setup(timestep=0.01)

Ecells = sim.Population(20, sim.HH_cond_exp())
Icells = sim.Population(10, sim.HH_cond_exp())
Net = Ecells + Icells
Esyn = sim.StaticSynapse(weight=0.05, delay=1.0)
random = sim.FixedProbabilityConnector(p_connect=0.5)
Econnections = sim.Projection(Ecells, Net, random, Esyn,
                              receptor_type='excitatory')
Isyn = sim.StaticSynapse(weight=0.2, delay=1.0)
Iconnections = sim.Projection(Icells, Net, random, Isyn,
                              receptor_type='inhibitory')

noise = sim.NoisyCurrentSource(mean=0.5, stdev=0.2, start=50.0,
                               stop=450.0, dt=1.0)
Ecells.inject(noise)
Net.record('v')
sim.run(1000.0)

data = Net.get_data()

sim.end()

arr = data.segments[0].analogsignalarrays[0]

# Plot E-cells
ax1 = pl.subplot(211)
for k in range(20):
    pl.plot(arr.times, arr[:, k], linewidth=2)
    pl.hold(True)

pl.xlabel('Time (ms)', fontsize=16)
pl.ylabel('Membrane potential (mV)', fontsize=16)
pl.title('E-Cells', fontsize=16)
pl.tick_params(labelsize=16)


# Plot I-cells
ax2 = pl.subplot(212, sharex=ax1)
for k in range(20, 30):
    pl.plot(arr.times, arr[:, k], linewidth=2)
    pl.hold(True)

pl.xlabel('Time (ms)', fontsize=16)
pl.ylabel('Membrane potential (mV)', fontsize=16)
pl.title('I-Cells', fontsize=16)
pl.tick_params(labelsize=16)
pl.show()
