# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 16:37:23 2013

@author: hari
"""
import numpy as np
from scipy import io, linalg
import pylab as pl

def plotDepthResults(subjlist, what,numCondsToPlot,summary = False, max = False,
                     ch_sel = 30):
        
    froot = '/home/hari/Documents/PythonCodes/research/DepthResults/'   
    #froot = '/home/hari/Documents/DepthResults/'   
    
    nsubjs = len(subjlist)
    condlist = [[1,7],[2,8],[3,9],[4,10]]
    condstemlist = ['_0dB','_m4dB','_m8dB','_m12dB']
    nconds = len(condlist)
    whatever_all = np.zeros((nsubjs,nconds))
    
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
                        whatever[condind] = np.max(dat[what][:,f_ind])
                    else:                   
                        f_pca = (np.logical_and(f > 80, f < 115)).squeeze()
                        C = np.cov(dat[what][:,f_pca])
                        lambdas, wts = linalg.eigh(C)
                        w = wts[:,-1]/(wts[:,-1]).sum()
                        whatever[condind] = np.dot(w,dat[what])[f_ind]
                else:
                    whatever[condind] = dat[what][ch_sel,f_ind]
                
        whatever_all[k,:] = whatever
        pl.plot(10*np.log10(whatever[0:numCondsToPlot]/1e-12),'o-',
                linewidth = 2)    
        pl.ylabel(what + ' (dB re: 1 sq.micro.V)', fontsize = 20)
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
    


# I36 was recorded with wrong sampling rate, I16 is missing
subjlist = ['I01','I02','I03','I05', 'I06','I07','I08','I09','I11','I13','I14',
            'I15','I17_redo','I18','I19','I20','I25','I26','I27','I28','I29',
            'I30','I32','I33','I34','I35','I37','I39','I16','I36', 'I41']
what = 'S'
numCondsToPlot = 4
S = plotDepthResults(subjlist,what,numCondsToPlot)
            