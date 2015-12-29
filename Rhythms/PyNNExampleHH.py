# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 14:36:30 2015

@author: Hari
"""
import pyNN.neuron as sim
import pylab as pl
from quantities import nA

sim.setup(timestep=0.01)

cell = sim.Population(1, sim.HH_cond_exp())
step_current = sim.DCSource(start=20.0, stop=150.0)
step_current.inject_into(cell)

cell.record('v')

for amp in [-0.2, 0.1]:
    step_current.amplitude = amp
    sim.run(200.0)
    sim.reset(annotations={"amplitude": amp*nA})

data = cell.get_data()

sim.end()

for segment in data.segments:
    vm = segment.analogsignalarrays[0]
    pl.plot(vm.times, vm, linewidth=2,
            label=str(segment.annotations["amplitude"]))
pl.legend(loc=0, fontsize=16)
pl.xlabel("Time (%s)" % vm.times.units._dimensionality, fontsize=16)
pl.ylabel("Membrane potential (%s)" % vm.units._dimensionality,
          fontsize=16)
pl.tick_params(labelsize=16)
pl.show()
