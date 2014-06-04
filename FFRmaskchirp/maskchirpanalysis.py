from anlffr.helper import biosemi2mne as bs
import mne
import numpy as np
from anlffr import spectral
from scipy import io
import os
import fnmatch

# Adding Files and locations
froot = '/home/hari/Documents/MATLAB/MaskChirp/'

# List of files stems, each will be appended by run number 
# Could be different for edf, fiff, eve etc.
# Use list [] and enumerate over if filenames are weird


# subjlist =['I01','I02','I03','I06','I07','I08','I09','I11','I13','I14','I15',
# 'I17_redo','I18','I19','I20','I25','I26','I27','I28','I29','I30','I37',
# 'I16','I32','I33','I34','I35','I39','I05','I36']

subjlist = ['I03']          
            
for subj in subjlist:
    
    fpath = froot + subj + '/'
    
    # These are so that the generated files are organized better
    respath = fpath + 'RES/'
    
    # Saving all the results (only) in a separate directory
    resonly_backup = '/home/hari/Documents/MaskChirpResults/' + subj + '/'
    
    #condlist = [[1,7],[2,8],[3,9],[4,10],[5,11]]
    condlist = [[6,12]]
    #condstemlist = ['_mask1','_mask2','_mask3','_mask4','_mask5']
    condstemlist = ['_chirp']
    
    for  condind, cond in enumerate(condlist):
        condstem = condstemlist[condind]
        print 'Running Subject', subj, 'Condition', condind
        
        save_raw_name = subj + condstem + '_alltrial.mat'
        
        if os.path.isfile(respath + save_raw_name):
            print 'Epoched data is already available on disk!'
            print 'Loading data from:', respath+save_raw_name
            x = io.loadmat(respath+save_raw_name)['x']
        else:
            bdfs = fnmatch.filter(os.listdir(fpath),subj+'*.bdf')
            print 'No pre-epoched data found, looking for BDF files'
            print 'Viola!', len(bdfs),  'files found!'
            
            for k, edfname in enumerate(bdfs):
                # Load data and read event channel
                (raw,eves) = bs.importbdf(fpath + edfname, nchans = 35,
                                          refchans = ['EXG1','EXG2'])
                
                raw.info['bads']+= ['A14','A25']
                # Filter the data
                raw.filter(l_freq = 70, h_freq = 1500, picks = np.arange(0,32,1))
                
                #raw.apply_proj()
                fs = raw.info['sfreq']
                
                # Here events 1 and 7 represent a particular stimulus in each polarity
                selectedEve = dict(up = cond[0], down = cond[1])
                
                # Epoching events of type 1 and 7
                epochs = mne.Epochs(raw,eves,selectedEve,tmin = -0.05, proj = False,
                                    tmax = 1.25, baseline = (-0.05, 0),
                                    reject = dict(eeg=150e-6))
                # Combining both polarities so I can get envelope related FFR responses
                epochs = mne.epochs.combine_event_ids(epochs,['up','down'],
                                                      dict(all= 101))
                # Getting the epoched data out, this step will also perform rejection
                xtemp = epochs.get_data()
                
                # Reshaping to the format needed by spectral.mtcpca() and calling it
                if(xtemp.shape[0] > 0):                 
                    xtemp = xtemp.transpose((1,0,2))
                    xtemp = xtemp[0:32,:,:]
                    
                    if(k==0):
                        x = xtemp
                    else:
                        x = np.concatenate((x,xtemp),axis = 1)
                else:
                    continue
                
    
    
        #params = dict(Fs=fs,fpass=[5,2000],tapers=[15, 29],Npairs = 2000,
        #              itc = 1)
        nPerDraw = 400
        nDraws = 100
        fs = 8096
        params = dict(Fs=fs,fpass=[5,2000],tapers=[15, 29],Npairs = 2000,
                    itc = 1)
        
        #        print 'Running Pairwise Spectrum Estimation'
        #       (pS,f) = spectral.mtpspec(x, params, verbose = 'DEBUG')
           
        print 'Running Raw Spectrum Estimation'
        (Sraw,f) = spectral.mtspecraw(x, params, verbose = True)
           
        print 'Running Mean Spectrum Estimation'
        (S,N,f) = spectral.mtspec(x,params,verbose = True)
           
        print 'Running CPCA PLV Estimation'
        (cplv,f)  = spectral.mtcpca(x, params, verbose = True)
        
        print 'Running CPCA PLV Estimation'
        (plv,f)  = spectral.mtplv(x, params, verbose = True)
        
        print 'Running CPCA PLV Estimation'
        (cpow,f)  = spectral.mtcspec(x, params, verbose = True)
        
        
        
        
        # Saving Results
        res = dict(cpow = cpow,plv = plv, cplv = cplv,
                   Sraw = Sraw, f = f, S = S, N = N)
        
        save_name = subj + condstem + '_results.mat'
        
        
        if (not os.path.isdir(respath)):
            os.mkdir(respath)
        if (not os.path.isdir(resonly_backup)):
            os.mkdir(resonly_backup)    
        io.savemat(respath + save_name,res)
        io.savemat(resonly_backup + save_name,res)
        
        if not os.path.isfile(respath + save_raw_name):
            io.savemat(respath + save_raw_name,dict(x = x,subj = subj))
        
'''
import pylab as pl
dat = io.loadmat(respath + subj + '_mask1_results.mat')
m0 = dat['cplv']
dat = io.loadmat(respath + subj + '_mask2_results.mat')
m = dat['cplv']
dat = io.loadmat(respath + subj + '_mask3_results.mat')
m_hi1 = dat['cplv']
dat = io.loadmat(respath + subj + '_mask4_results.mat')
m_lo = dat['cplv']
dat = io.loadmat(respath + subj + '_mask5_results.mat')
m_hi2 = dat['cplv']
pl.plot(f,m0)
pl.hold(True)
pl.plot(f,m)
pl.plot(f,m_hi1)
pl.plot(f,m_lo)
pl.plot(f,m_hi2)
pl.legend(('Quiet','On Freq','High 1','Low', 'High 2'))       
pl.show()
'''
import pylab as pl
pl.plot(f,cplv)
pl.show()
