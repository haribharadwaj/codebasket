# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from scipy import io
import numpy as np
import pylab as pl


dat = io.loadmat('ERPsummary_zscore.mat')
t = dat['t'].flatten()
c7 = dat['c7']
c14 = dat['c14']
c20 = dat['c20']

start, stop = 0.1, 0.4

s7 = c7[:, np.logical_and(t > start, t < stop)].mean(axis=1)
s14 = c14[:, np.logical_and(t > start, t < stop)].mean(axis=1)
s20 = c20[:, np.logical_and(t > start, t < stop)].mean(axis=1)

x = np.asarray([5.5, 11.5, 17.5])
y = np.asarray([s7[:26].mean(), s14[:26].mean(), s20[:26].mean()])
yerr = np.asarray([s7[:26].std() / (26 ** 0.5), s14[:26].std() / (26 ** 0.5),
                   s20[:26].std() / (26 ** 0.5)])
pl.errorbar(x, y, yerr,
            fmt='ok-', elinewidth=2)

x = x + 1
y = np.asarray([s7[26:].mean(), s14[26:].mean(), s20[26:].mean()])
yerr = np.asarray([s7[26:].std() / (21 ** 0.5), s14[26:].std() / (21 ** 0.5),
                   s20[26:].std() / (21 ** 0.5)])
pl.errorbar(x, y, yerr,
            fmt='sr--', elinewidth=2)
pl.xlabel('Number of Coherent Tones', fontsize=16)
pl.ylabel('Evoked Response (normalized)', fontsize=16)
pl.xticks((6, 12, 18))
ax = pl.gca()
ax.tick_params(labelsize=14)
pl.legend(('TD', 'ASD'), loc=0)

csvfile = open('ERPforR.csv', 'w')
csvfile.write('subj, group, erp, cond\n')
for k in range(26):
    csvfile.write(str(k+1) + ', TD,' + str(s7[k]) + ', low\n')
    csvfile.write(str(k+1) + ', TD,' + str(s14[k]) + ', mid\n')
    csvfile.write(str(k+1) + ', TD,' + str(s20[k]) + ', high\n')

for k in range(26, 47):
    csvfile.write(str(k+1) + ', ASD,' + str(s7[k]) + ', low\n')
    csvfile.write(str(k+1) + ', ASD,' + str(s14[k]) + ', mid\n')
    csvfile.write(str(k+1) + ', ASD,' + str(s20[k]) + ', high\n')

csvfile.close()