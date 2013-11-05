# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 22:17:31 2013

Script to analyze BW data and get Qerb
@author: Hari Bharadwaj
"""
from scipy import io
from glob import glob
import pylab as pl
import numpy as np


rootdir = '/home/hari/Documents/PythonCodes/research/BW/'
subj = 'I25'

flist = glob(rootdir + subj + '/*.mat')
bwlist = []
threshlist = []
asymmlist = []
for k, fname in enumerate(flist):
    dat = io.loadmat(fname)   
    if(dat.has_key('bw')):
        bw = dat['bw'].ravel()[0]
        asymm = dat['asymm'].ravel()[0]
        thresh = dat['thresh'].ravel()[0] - dat['soundLevel'].ravel()[0]
        
        bwlist = bwlist + [bw,]
        asymmlist = asymmlist + [asymm,]
        threshlist = threshlist + [thresh,]
        if(asymm == 0):
            pl.plot(bw,thresh,'ok',linewidth = 3, ms = 10)
            pl.hold(True)
        elif(asymm == 1):
            pl.plot(bw,thresh,'<b',linewidth = 3, ms = 10)
            pl.hold(True)
        elif(asymm == 2):
            pl.plot(bw,thresh, '>r',linewidth = 3, ms = 10)
            pl.hold(True)
        else:
            print 'Unknown notch shape!'
            
bwlist = np.asarray(bwlist)
threshlist = np.asarray(threshlist)
asymmlist = np.asarray(asymmlist)    
pl.xlabel('Notch Width',fontsize = 20)
pl.ylabel('Threshold (dB SL)',fontsize = 20)
pl.xlim((-0.02, 0.22))
pl.show()


    
    
    