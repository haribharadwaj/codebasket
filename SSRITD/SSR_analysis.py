from anlffr.helper import biosemi2mne as bs
import mne
import numpy as np
from anlffr import spectral
from scipy import io
import os
import fnmatch

# Adding Files and locations
froot = '/home/hari/Documents/MATLAB/SSRITD/EEG/'

# List of files stems, each will be appended by run number 
# Could be different for edf, fiff, eve etc.
# Use list [] and enumerate over if filenames are weird


subjlist = ['I02','I03','I05','I06','I08','I09','I11','I13','I14','I15','I17',
            'I18','I19','I20','I25','I26','I29','I30','I33','I35','I36','I37'] 
            
for subj in subjlist:
    
    fpath = froot + subj + '/'
    
    # These are so that the generated files are organized better
    respath = fpath + 'RES/'
    
    # Saving all the results (only) in a separate directory
    resonly_backup = '/home/hari/Documents/SSRresults/' + subj + '/'
    
    condlist = [[8,9,10,11,12]]
    condstemlist = ['_allITDs']
    
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
                (raw,eves) = bs.importbdf(fpath + edfname)
                
                # Filter the data
                raw.filter(l_freq = 10, h_freq = 200, picks = np.arange(0,32,1))
                
                conds = {'50us':8,'100us':9,'200us':10,'400us':11,'800us':12}
                selectedEves = ['8','9','10','11','12'] # MNE bug??
               
                # Epoching events of type 1 and 7
                epochs = mne.Epochs(raw,eves,cond,tmin = 1.2, proj = False,
                                    tmax = 2.0, baseline = (1.2, 2.0),
                                    reject = dict(eeg=150e-6))
                # Combining both polarities so I can get envelope related FFR responses
                epochs = mne.epochs.combine_event_ids(epochs,selectedEves,
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
                
    
    
        params = dict(Fs=2048.0,fpass=[5,600],tapers=[1, 1],pad=1,Npairs = 2000,
                      itc = 1)
        
        
        print 'Running Pairwise Spectrum Estimation'
        (pS,f) = spectral.mtpspec(x, params, verbose = 'DEBUG')
        
        print 'Running Raw Spectrum Estimation'
        (Sraw,f) = spectral.mtspecraw(x, params, verbose = True)
        
        print 'Running Mean Spectrum Estimation'
        (S,N,f) = spectral.mtspec(x,params,verbose = True)
        
        print 'Running CPCA PLV Estimation'
        (cplv,f)  = spectral.mtcpca(x, params, verbose = True)
        
        
    
        
    
        # Saving Results
        res = dict(pS = pS,cplv = cplv,Sraw = Sraw, f = f, S = S, N = N)
        save_name = subj + condstem + '.mat'
        
        
        if (not os.path.isdir(respath)):
            os.mkdir(respath)
        if (not os.path.isdir(resonly_backup)):
            os.mkdir(resonly_backup)    
        io.savemat(respath + save_name,res)
        io.savemat(resonly_backup + save_name,res)
        
        if not os.path.isfile(respath + save_raw_name):
            io.savemat(respath + save_raw_name,dict(x = x,subj = subj))
           

