# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 16:37:23 2013

@author: hari
"""
import numpy as np
from scipy import io, linalg
import pylab as pl

def plotSSRResults(subjlist, what,summary = False, max = False,
                     ch_sel = 30):
        
    froot = '/home/hari/Documents/PythonCodes/research/SSRITD/EEG/'   
    
    nsubjs = len(subjlist)
    
    whatever_all = np.zeros((nsubjs))
    
    for k, subj in enumerate(subjlist):
        fpath = froot + subj + '/'
        load_name = fpath + subj + '_allITDs.mat'
        dat = io.loadmat(load_name)
        f = dat['f']
        f_ind = np.argmin(abs(f - 40))
        
        if(what == 'cplv'):
            summary = False
            whatever = dat[what][f_ind]
        
        else:
                    
            if(summary):
                if(max):
                    whatever = np.max(dat[what][:,f_ind])
                else:                   
                    f_pca = (np.logical_and(f > 35, f < 45)).squeeze()
                    C = np.cov(dat[what][:,f_pca])
                    lambdas, wts = linalg.eigh(C)
                    w = wts[:,-1]/(wts[:,-1]).sum()
                    whatever = np.dot(w,dat[what])[f_ind]
            else:
                whatever = dat[what][ch_sel,f_ind]
            
        whatever_all[k] = whatever
        pl.plot(20*np.log10(whatever/1e-6),'o-',
                linewidth = 2)    
        pl.ylabel(what + ' (dB re: 1 uV)', fontsize = 20)
        pl.xlabel('Modulation Depth',fontsize = 20)
        pl.hold(True)
        ax = pl.gca()
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20)    
        
    pl.show()
    pl.legend(subjlist)        
    return  whatever_all
    


subjlist = ['I02','I03','I05','I06','I08','I09','I11','I13','I14','I15','I17',
            'I18','I19','I20','I25','I26','I29','I30','I33','I35','I36','I37'] 
  
what = 'cplv'
cplv = plotSSRResults(subjlist,what)
            