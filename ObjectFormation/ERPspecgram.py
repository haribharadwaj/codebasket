# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from scipy import io
from anlffr import tfr
import numpy as np


dat = io.loadmat('ERPsummary_zscore.mat')
t = dat['t'].flatten()
c7 = dat['c7']
c14 = dat['c14']
c20 = dat['c20']

f = np.arange(5, 200)
powertd, itc, times = tfr.tfr_multitaper(c20[:26, None, :], 1000., f)
times = times + t[0]
powertd_adj = tfr.rescale(powertd, times, (-0.2, 0.), mode='logratio')
tfr.plot_tfr(powertd_adj, times, f)


powerasd, itc, times = tfr.tfr_multitaper(c20[26:, None, :], 1000., f)
times = times + t[0]
powerasd_adj = tfr.rescale(powerasd, times, (-0.2, 0.), mode='logratio')
tfr.plot_tfr(powerasd_adj, times, f)
