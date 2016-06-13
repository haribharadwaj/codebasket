# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 21:16:41 2016

@author: Hari
"""

import numpy as np


class singleunit:
    def __init__(self, b=-0.01, s=0, theta=0., label=None):
        self.theta = theta
        self.b = b
        self.s = s
        self.label = label

    def update(self, b, s):
        self.s = s
        self.b = b

    def step(self, dt=0.01):
        I = self.b + self.s
        dtheta = ((1 - np.cos(self.theta)) + I * (1 + np.cos(self.theta))) * dt
        self.theta = self.theta + dtheta
        if self.theta > np.pi:
            self.theta = self.theta - 2 * np.pi


if __name__ == "__main__":
    dt = 0.01
    Ntime = 10000
    u = np.ones(Ntime) * 0.25
    t = np.arange(Ntime) * dt
    c = singleunit()
    th = np.zeros(t.shape)
    for k, su in enumerate(u):
        c.update(c.b, su)
        c.step()
        th[k] = c.theta
    import pylab as pl
    pl.plot(t, th)
    pl.xlabel('Time (ms)')
    pl.ylabel('Neuron State (theta)')
