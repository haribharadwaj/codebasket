# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 12:35:30 2014

@author: hari
"""
import numpy as np
from scipy import io
import pylab as pl

subjlist = ['I33', 'I41', 'I08', 'I02', 'I11', 'I13', 'I52', 'I14']

# Adding Files and locations
froot = '/home/hari/Documents/PythonCodes/ModelData/'

pl.figure(num=1)
pl.figure(num=2)
Pe = 0
P2e = 0
PLVe = 0
PLV2e = 0
Phe = 0
Ph2e = 0
te = 0
t2e = 0

for subj in subjlist:
    respath = froot + subj + '/RES/'
    # Magnitude plots
    condlist = np.arange(1, 16)
    f = (condlist-1)*30 + 100
    condstemlist = np.array(map(str, f))
    plv = np.zeros(condlist.shape)
    P = np.zeros(condlist.shape)
    for k, cond in enumerate(condstemlist):
        fname = respath + subj + '_' + cond + 'Hz_results.mat'
        f0 = f[k]
        dat = io.loadmat(fname)
        f_fft = dat['f'].squeeze()
        P[k] = 10 * (np.log10(dat['cpow'].squeeze()
                     [np.argmin(abs(f_fft-f0))]/1e-12))
        plv[k] = dat['cplv'].squeeze()[np.argmin(abs(f_fft-f0))]

    pl.figure(num=1)
    ax1 = pl.subplot(2, 1, 1)
    pl.hold(True)
    Pe = Pe + P*(1.0/len(subjlist))
    P2e = P2e + (P**2)*(1.0/len(subjlist))
    pl.plot(f, P, linewidth=2)
    pl.ylabel('Power (dB re: 1uV)', fontsize=16)
    pl.subplot(2, 1, 2, sharex=ax1)
    pl.hold(True)
    PLVe = PLVe + plv*(1.0/len(subjlist))
    PLV2e = PLV2e + (plv**2)*(1.0/len(subjlist))
    pl.plot(f, plv, linewidth=2)
    pl.ylabel('PLV (squared)', fontsize=16)
    pl.xlabel('Modulation frequency (Hz)', fontsize=16)

    # Phase plots
    fname = respath + subj + '_all_phase.mat'
    datph = io.loadmat(fname)
    ph = datph['Ph_f0'].squeeze()
    t0 = -0.025
    shift = 2*np.pi*f*t0
    ph = ph - shift
    # phclean = np.median(unwrap(ph, axis=1), axis=0) / (2 * np.pi)
    phclean = np.unwrap(ph[24, :]) / (2 * np.pi)
    phclean = phclean - phclean.mean()
    f = datph['f0_list'].squeeze()
    offset = -1.6
    t_grp = -1e3 * np.diff(phclean) / (np.diff(f).mean()) + offset
    f_grp = f[1:] - (np.diff(f).mean())/2

    pl.figure(num=2)
    ax2 = pl.subplot(2, 1, 1)
    pl.hold(True)
    Phe = Phe + phclean*(1.0/len(subjlist))
    Ph2e = Ph2e + (phclean**2)*(1.0/len(subjlist))
    pl.plot(f, phclean, linewidth=2)
    pl.ylabel('Unwrapped phase (cycles)', fontsize=16)
    pl.subplot(2, 1, 2, sharex=ax2)
    te = te + t_grp*(1.0/len(subjlist))
    t2e = t2e + (t_grp**2)*(1.0/len(subjlist))
    pl.plot(f_grp, t_grp, linewidth=2)
    pl.hold(True)
    pl.xlabel('Modulation frequency (Hz)', fontsize=16)
    pl.ylabel('Group Delay (ms)', fontsize=16)
    pl.xlim(np.min(f), np.max(f))
    pl.legend(subjlist, loc=1)

Perr = ((P2e - Pe**2) / len(subjlist))**0.5
PLVerr = ((PLV2e - PLVe**2) / len(subjlist))**0.5
Pherr = ((Ph2e - Phe**2) / len(subjlist))**0.5
terr = ((t2e - te**2) / len(subjlist))**0.5

pl.show()

pl.figure(num=3)
ax3 = pl.subplot(2, 2, 1)
pl.plot(f, Pe, linewidth=2)
pl.hold(True)
pl.plot(f, Pe + Perr, 'r--')
pl.plot(f, Pe - Perr, 'r--')
pl.ylabel('Power (dB re: 1uV)', fontsize=16)
pl.subplot(2, 2, 3, sharex=ax3)
pl.plot(f, PLVe, linewidth=2)
pl.hold(True)
pl.plot(f, PLVe + PLVerr, 'r--')
pl.plot(f, PLVe - PLVerr, 'r--')
pl.ylabel('PLV (squared)', fontsize=16)
pl.xlabel('Modulation frequency (Hz)', fontsize=16)
pl.subplot(2, 2, 2)
pl.plot(f, Phe, linewidth=2)
pl.hold(True)
pl.plot(f, Phe + Pherr, 'r--')
pl.plot(f, Phe - Pherr, 'r--')
pl.ylabel('Unwrapped phase (cycles)', fontsize=16)
pl.subplot(2, 2, 4, sharex=ax3)
te_new = -1e3 * np.diff(Phe) / (np.diff(f).mean()) + offset
pl.plot(f_grp, te_new, linewidth=2)
pl.hold(True)
pl.plot(f_grp, te + terr, 'r--')
pl.plot(f_grp, te - terr, 'r--')
pl.ylim((-5, 15.))
pl.xlabel('Modulation frequency (Hz)', fontsize=16)
pl.ylabel('Group Delay (ms)', fontsize=16)
pl.show()
