# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 12:35:30 2014

@author: hari
"""
import numpy as np
from scipy import io
import pylab as pl

subj = 'I14'

# Adding Files and locations
froot = '/home/hari/Documents/PythonCodes/ModelData/'

respath = froot + subj + '/RES/'
# Magnitude plots
condlist = np.arange(1, 13)  # Only plotting upto 430 Hz
f = (condlist-1)*30 + 100
condstemlist = np.array(map(str, f))

for k, cond in enumerate(condstemlist):
    print k
    fname = respath + subj + '_' + cond + 'Hz_results.mat'
    f0 = f[k]
    dat = io.loadmat(fname)
    f_fft = dat['f'].squeeze()
    P = 10 * (np.log10(dat['cpow'].squeeze() / 1e-12))
    plv = dat['cplv'].squeeze()

    pl.figure(num=1)
    ax1 = pl.subplot(4, 3, k+1)
    pl.hold(True)
    pl.plot(f_fft, P, linewidth=2)
    pl.ylabel('Power (dB re: 1uV)', fontsize=16)

    pl.figure(num=2)
    pl.hold(True)
    ax1 = pl.subplot(4, 3, k+1)
    pl.plot(f_fft, plv, linewidth=2)
    pl.ylabel('PLV (squared)', fontsize=16)
    pl.xlabel('Frequency (Hz)', fontsize=16)
    pl.title('fm = ' + str(f0) + ' Hz')
pl.show()
