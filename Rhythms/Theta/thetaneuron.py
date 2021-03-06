# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 21:16:41 2016

@author: Hari
"""

import numpy as np


class singleunit:
    """Single theta neuron

    Parameters
    ----------
    theta : double
        State variable.
    b : double
        Bias parameter (determines stable state in the absence of input).
    s : double
        Instantaneous synaptic input.
    label : str, int, or None.
        Label assigned to identify the unit.

    Attributes
    ----------
    theta : double
        State variable.
    b : double
        Bias parameter (determines stable state in the absence of input).
    s : double
        Instantaneous synaptic input.
    label : str, int, or None.
        Label assigned to identify the unit.

    Notes
    -----
    The theta neuron is a canonical model neuron exhibiting a saddle-node
    bifurcation on an invariant cycle. It can be obtained from a quadratic
    integrate-and-fire neuron through a change of variables and from a Hodgkin-
    Huxley type-I neuron through dimension reduction approximations.

    [1] Ermentrout, G. B., & Kopell, N. (1986). Parabolic bursting in an
    excitable system coupled with a slow oscillation. SIAM Journal on Applied
    Mathematics, 46(2), 233-253.

    [2] Vierling-Claassen, D., Siekmeier, P., Stufflebeam, S., & Kopell, N.
    (2008). Modeling GABA alterations in schizophrenia: a link between impaired
    inhibition and altered gamma and beta range auditory entrainment. Journal
    of neurophysiology, 99(5), 2656-2671.

    [3] Börgers, C., & Kopell, N. J. (2008). Gamma oscillations and stimulus
    selection. Neural computation, 20(2), 383-414.

    """
    def __init__(self, b=-0.01, s=0., theta=0., label=None):
        self.theta = theta
        self.b = b
        self.s = s
        self.label = label

    def update(self, s, b=None):
        """Update inputs and/or the bias parameter of the theta neuron

        Parameters
        ----------
        s : float
            Instantaneous synaptic drive
        b : float
            Bias parameter of the theta neuron

        Returns
        -------
        None : The attributes are modified in place

        """
        self.s = s
        if b is not None:
            self.b = b

    def step(self, dt=0.01):
        """Simulate a single forward Euler step for the theta neuron state

        Parameters
        ----------
        dt : double
            The simulation time step in milliseconds

        Returns
        -------
        None : The state is modified in place

        """
        I = self.b + self.s
        dtheta = ((1 - np.cos(self.theta)) + I * (1 + np.cos(self.theta))) * dt
        self.theta = self.theta + dtheta
        if self.theta > np.pi:
            self.theta = self.theta - 2 * np.pi


class synapse:
    def __init__(self, tr=0.1, td=8., m=0., eta=5.):
        self.tr = tr
        self.td = td
        self.m = m
        self.eta = eta

    def step(self, pre, dt=0.01):
        dm = (-self.m / self.td +
              np.exp(-self.eta * (1 + np.cos(pre.theta))) *
              (1 - self.m) / self.tr) * dt
        self.m = self.m + dm
        if self.m > 1.:
            self.m = 1.
        if self.m < 0.:
            self.m = 0.


class network:
    """Theta neuron network

    Parameters
    ----------
    units : list of singleunits
        The collection of neurons that form the network.
    exOrInh : list of int (same length as units)
        Whether each cell is excitatory (+1) or inhibitory (-1)
    cMat : boolean numpy.2darray
        Connectivity matrix (only wiring, i.e., connected or not).
    gee : double
        Excitatory on to Excitatory (E-E) synaptic strength (homogenous).
    gei : double
        E-I synaptic strength (homogeneous).
    gie : double
        I-E synaptic strength (homogenous).
    gii : double
        I-I synaptic strength (homogeneous).

    Attributes
    ----------
    units : list of singleunits or int
        The collection of neurons that form the network, or number of neurons.
        If latter, singleunits with default attributes and integer labels are
        instantiated and placed in a list which is then assigned to the units
        attribute of the network.
    exOrInh : list of boolean (same length as units)
        Whether each cell is excitatory (False) or inhibitory (True)
    cMat : boolean numpy.2darray
        Connectivity matrix (only wiring, i.e., connected or not).
    gee : double
        Excitatory on to Excitatory (E-E) synaptic strength (homogenous).
    gei : double
        E-I synaptic strength (homogeneous).
    gie : double
        I-E synaptic strength (homogenous).
    gii : double
        I-I synaptic strength (homogeneous).

    Notes
    -----
    This is a simple homogenous E-I network of theta neurons built based on the
    simplified model of Vierling-Claassen et al., (2008).

    [1] Vierling-Claassen, D., Siekmeier, P., Stufflebeam, S., & Kopell, N.
    (2008). Modeling GABA alterations in schizophrenia: a link between impaired
    inhibition and altered gamma and beta range auditory entrainment. Journal
    of neurophysiology, 99(5), 2656-2671.

    """
    def __init__(self, units, exOrInh, cMat, gee=0.015, gei=0.025, gie=0.015,
                 gii=0.02):
        if type(units) is list:
            self.units = units
        else:
            self.units = []
            for k in range(units):
                self.units += [singleunit(label=k)]
        self.cMat = cMat
        self.exOrInh = exOrInh
        self.gee = gee
        self.gei = gei
        self.gie = gie
        self.gii = gii
        self._syn = []
        self._Ncells = len(self.units)
        for j in range(self._Ncells):
            temp = []
            for k in range(self._Ncells):
                temp += [synapse(td=8.) if self.exOrInh[j] == -1.
                         else synapse(td=2.)]
            self._syn += [temp, ]

    def step(self, dt=0.01):
        for k in range(self._Ncells):  # Loop over post-synaptic cells
            for j in range(self._Ncells):  # Loop over pre-synaptic cells
                self._syn[j][k].step(self.units[j], dt=dt)
                # Add synapting input from other cells in the network
                if self.exOrInh[j] == 1.:
                    if self.exOrInh[k] == 1.:
                        g = self.gee
                    else:
                        g = self.gei
                else:
                    if self.exOrInh[k] == 1.:
                        g = self.gie
                    else:
                        g = self.gii
                self.units[k].s += (self.exOrInh[j] * self.cMat[j, k] *
                                    self._syn[j][k].m * g)
        for k in range(self._Ncells):  # Loop again separately
            self.units[k].step(dt=dt)
            self.units[k].s = 0.


if __name__ == "__main__":
    dt = 0.2
    tmax = 1000.
    t = np.arange(0., tmax, dt)
    # Initialize network
    Ncells = 30
    Nex = np.int(Ncells * 2/3)
    Ninh = Ncells - Nex
    exOrInh = [1] * Nex + [-1] * Ninh

    case = 'sine'

    ur = np.random.randn(t.shape[0]) * 0.15 + 0.2
    ur[t < 50.] = 0.
    ur[t > 800.] = 0.

    fs = 25.
    us = np.zeros(t.shape[0])
    ts_sine = np.arange(20., tmax, 1000./fs)
    pulse = 0.5 * (np.exp(-t/0.1) - np.exp(-t/2.))/(-1.9)
    for sp_sine in ts_sine:
        us[np.argmin(np.abs(t - sp_sine))] = 1.
    us = np.convolve(us, pulse)[:t.shape[0]]
    us[t < 100.] = 0.
    # us[t > 800.] = 0.

    # Noise input for each cell independently
    inc = 30.
    Nspikes = np.int(np.ceil(tmax / inc))
    ts = np.cumsum(np.random.exponential(scale=inc, size=(Ncells, Nspikes)),
                   axis=1)
    noise = np.zeros((Ncells, t.shape[0]))
    for cell in range(Ncells):
        for sp in range(Nspikes):
            noise[cell, np.argmin(np.abs(t - ts[cell, sp]))] = 1.
    for cell in range(Ncells):
        noise[cell, :] = np.convolve(noise[cell, :], pulse)[:t.shape[0]]

    driverin = singleunit()
    driverout = synapse(td=2.)
    inputlist = dict(random=ur, sine=us, spont=np.zeros(t.shape[0]))
    u = inputlist[case]

    # Connectivity matrix and initialize network
    cMat = np.random.rand(Ncells, Ncells) > 0.
    for k in range(Ncells):
        cMat[k, k] = 0
    N = network(Ncells, exOrInh, cMat)  # gei=0.1, gii=0.05, gie=0.01, gee=0.)

    gde = 0.3 / pulse.max()
    gdi = 0.08 / pulse.max()
    gne = 0.7
    gni = 0.7

    # Choose cellw to input noise
    inputCell = range(Ncells)

    th = np.zeros((Ncells, t.shape[0]))
    curr = np.zeros((Ncells, t.shape[0]))

    # Run simulation
    for k, uk in enumerate(u):
        if np.mod(k, t.shape[0]/20) == 0:
            print 'time = %f / %f' % (t[k], t[-1])

        for i in range(Nex):
            if i in inputCell:
                N.units[i].update(s=uk * gde + noise[i, k] * gne)
            else:
                N.units[i].update(s=noise[i, k])
        for i in range(Ninh):
            if i+Nex in inputCell:
                N.units[i + Nex].update(s=(uk * gdi +
                                           noise[i + Nex, k] * gni))
            else:
                N.units[i + Nex].update(s=noise[i + Nex, k])

        N.step(dt=dt)
        for plotCell in range(Ncells):
            th[plotCell, k] = N.units[plotCell].theta
            for j in range(Ncells):
                if N.exOrInh[j] == 1.:
                    if N.exOrInh[plotCell] == 1.:
                        g = N.gee
                    else:
                        g = N.gei
                else:
                    if N.exOrInh[plotCell] == 1.:
                        g = N.gie
                    else:
                        g = N.gii
                curr[plotCell, k] += (N._syn[j][plotCell].m * N.exOrInh[j] *
                                      g * N.cMat[j, plotCell])

    # Plot results
    import pylab as pl
    ax1 = pl.subplot(4, 1, 1)
    pl.plot(t, th[:20].T, 'b')
    pl.hold(True)
    pl.plot(t, th[20:].T, 'r')
    pl.ylabel('Neuron State (theta)')
    ax2 = pl.subplot(4, 1, 2, sharex=ax1)
    pl.plot(t, curr[:20].T, 'b')
    pl.hold(True)
    pl.plot(t, curr[20:].T, 'r')
    pl.ylabel('Input Current')
    ax3 = pl.subplot(4, 1, 3, sharex=ax1)
    pl.plot(t, u)
    pl.ylabel('Input Drive')
    pl.hold(True)
    ax4 = pl.subplot(4, 1, 4, sharex=ax1)
    MEG = curr[:20, ].sum(axis=0)
    pl.plot(t, MEG)
    pl.xlabel('Time (ms)')
    pl.show()
