# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 16:37:23 2013

@author: hari
"""
import numpy as np
from scipy import io
import pylab as pl

def plotDepthResults(subjlist, numCondsToPlot, harm = 1):
        
    #froot = '/home/hari/Documents/PythonCodes/research/DepthResults/'   
    froot = '/home/hari/Documents/DepthResults/'   
    
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
            load_name = fpath + subj + condstem + 'cpow.mat'
            dat = io.loadmat(load_name)
            f = dat['f']
            
            f_ind = np.argmin(abs(f - 100*harm))
            whatever[condind] = dat['cpow'][f_ind]
            
            
                
        whatever_all[k,:] = whatever
        pl.plot(10*np.log10(whatever[0:numCondsToPlot]/1e-12),'o-',
                linewidth = 2)    
        pl.ylabel('cPCA Power (dB re: 1 sq.micro.V)', fontsize = 20)
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
    

subjlist = ['I01','I02','I03','I05', 'I06','I07','I08','I09','I11','I13','I14',
            'I15','I16','I17_redo','I18','I19','I20','I25','I26','I27','I28',
            'I29','I30','I32','I33','I34','I35','I37','I39','I16','I36']
numConds = 4

numHarms = 7
harms = np.arange(1,numHarms+1)

numSubjs = len(subjlist)

cS = np.zeros((numSubjs,numConds,numHarms))
for k, harm in enumerate(harms):
    cS[:,:,k] = plotDepthResults(subjlist,numConds, harm)

respE = cS.sum(axis = 2)
