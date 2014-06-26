# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 16:37:23 2013

@author: hari
"""
import numpy as np
from scipy import io, linalg
import pylab as pl


def plotDepthResults(subjlist, what, numCondsToPlot, summary=False, max=False,
                     ch_sel=30):

    froot = '/home/hari/Documents/PythonCodes/research/DepthResults/'
    #froot = '/home/hari/Documents/DepthResults/'

    nsubjs = len(subjlist)
    condlist = [[1, 7], [2, 8], [3, 9]]
    condstemlist = ['_0dB', '_m4dB', '_m8dB']
    nconds = len(condlist)
    whatever_all = np.zeros((nsubjs, nconds))

    for k, subj in enumerate(subjlist):
        fpath = froot + subj + '/'

        whatever = np.zeros(len(condlist))

        for condind, cond in enumerate(condlist):
            condstem = condstemlist[condind]
            load_name = fpath + subj + condstem + '.mat'
            dat = io.loadmat(load_name)
            f = dat['f']

            f_ind = np.argmin(abs(f - 100))

            if(what == 'cplv'):
                summary = False
                whatever[condind] = dat[what][f_ind]

            else:

                if(summary):
                    if(max):
                        whatever[condind] = np.max(dat[what][:, f_ind])
                    else:
                        f_pca = (np.logical_and(f > 80, f < 115)).squeeze()
                        C = np.cov(dat[what][:, f_pca])
                        lambdas, wts = linalg.eigh(C)
                        w = wts[:, -1] / (wts[:, -1]).sum()
                        whatever[condind] = np.dot(w, dat[what])[f_ind]
                else:
                    whatever[condind] = dat[what][ch_sel, f_ind]

        whatever_all[k, :] = whatever
        pl.plot(20 * np.log10(whatever[0:numCondsToPlot] / 1e-6), 'o-',
                linewidth=2)
        pl.ylabel(what + ' (dB re: 1 sq.micro.V)', fontsize=20)
        pl.xlabel('Modulation Depth', fontsize=20)
        pl.hold(True)
        ax = pl.gca()
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20)

    pl.show()
    pl.legend(subjlist)
    return whatever_all


# I36 was recorded with wrong sampling rate, I16 is missing
top = ['I37', 'I15', 'I39', 'I26', 'I03', 'I13', 'I33', 'I19', 'I02',
       'I30', 'I41', 'I14', 'I08']

bottom = ['I09', 'I18', 'I25', 'I05', 'I29',
          'I20', 'I36', 'I06', 'I11', 'I07', 'I05']

what = 'S'
numCondsToPlot = 3
S_top = np.log10(plotDepthResults(top, what, numCondsToPlot)/1e-6) * 20.0
S_bottom = np.log10(plotDepthResults(bottom, what, numCondsToPlot)/1e-6) * 20.0

m = [0, -4.0, -8.0]
top_mu = S_top.mean(axis=0)
top_err = S_top.std(axis=0)/(len(top)**0.5)

bottom_mu = S_bottom.mean(axis=0)
bottom_err = S_bottom.std(axis=0)/(len(bottom)**0.5)

pl.close('all')
pl.figure()
pl.errorbar(m, top_mu, yerr=top_err, linewidth=3)
pl.hold(True)
pl.errorbar(m, bottom_mu, yerr=bottom_err, color='r', linewidth=3)
pl.plot(m, S_top.T, '--', linewidth=2, color=[0.5, 0.5, 0.5])
pl.plot(m, S_bottom.T, '--', linewidth=2, color=[0.5, 0.5, 0.5])
pl.xlabel('Modulation Depth (dB re: 100%)', fontsize=20)
pl.ylabel('EFR magnitude (dB re: 1uV)', fontsize=20)
pl.xlim((-8.5, 0.5))
pl.show()
