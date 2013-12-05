from anlffr.helper import biosemi2mne as bs
import mne
import numpy as np
from anlffr import spectral
import pylab as pl
from scipy import io

# Adding Files and locations
froot = '/home/hari/Documents/MATLAB/Depth_Clean/'

# List of files stems, each will be appended by run number 
# Could be different for edf, fiff, eve etc.
# Use list [] and enumerate over if filenames are weird

subj = 'I08'
namestem = subj + '_depth'

fpath = froot + subj + '/'

# These are so that the generated files are organized better
respath = fpath + 'RES/'

nruns = 5
condlist = [[1,7],[2,8],[3,9],[4,10]]
condstemlist = ['_0dB','_m4dB','_m8dB','_m12dB']

for  condind, cond in enumerate(condlist):
    condstem = condstemlist[condind]
    print 'Running Subject', subj, 'Condition', condind
    for k in np.arange(0,nruns):
        # Load data and read event channel
        edfname = fpath + namestem + '_0' + str(k+1) + '.bdf'
        (raw,eves) = bs.importbdf(edfname)
        
        # Filter the data
        raw.filter(l_freq = 70, h_freq = 1500, picks = np.arange(0,32,1))
        
        # Here events 1 and 7 represent a particular stimulus in each polarity
        selectedEve = dict(up = cond[0], down = cond[1])
        
        # Epoching events of type 1 and 7
        epochs = mne.Epochs(raw,eves,selectedEve,tmin = -0.05, proj = False,
                            tmax = 0.45, baseline = (-0.05, 0),
                            reject = dict(eeg=100e-6))
        # Combining both polarities so I can get envelope related FFR responses
        epochs = mne.epochs.combine_event_ids(epochs,['up','down'],
                                              dict(all= 101))
        # Getting the epoched data out, this step will also perform rejection
        xtemp = epochs.get_data()
        
        # Reshaping to the format needed by spectral.mtcpca() and calling it
        xtemp = xtemp.transpose((1,0,2))
        xtemp = xtemp[0:32,:,:]
        
        if(k==0):
            x = xtemp
        else:
            x = np.concatenate((x,xtemp),axis = 1)


    params = dict(Fs=4096,fpass=[5,600],tapers=[1, 1],pad=1,Npairs = 2000,
                  itc = 1)
    nPerDraw = 400
    nDraws = 100
    
    
    print 'Running Pairwise Spectrum Estimation'
    (pS,f) = spectral.mtpspec(x, params)
    
    print 'Running Raw Spectrum Estimation'
    (Sraw,f) = spectral.mtspecraw(x, params)
    
    print 'Running CPCA PLV Estimation'
    (cplv,f)  = spectral.mtcpca(x, params)

    # Plotting results
    pl.plot(f,pS[30,:],linewidth = 2)
    pl.hold(True)
    

    # Saving Results
    res = dict(pS = pS,selectedEve = selectedEve, cplv = cplv,
               Sraw = Sraw, f = f)
    save_name = subj + condstem + '.mat'
    save_raw_name = subj + condstem + '_alltrial.mat'
    io.savemat(respath + save_name,res)

    io.savemat(respath + save_raw_name,
               dict(x = x, selectedEve = selectedEve, subj = subj))

pl.ylabel('Pairwise Power', fontsize = 20)
pl.xlabel('Frequency (Hz)',fontsize = 20)
pl.show()
ax = pl.gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(20)
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(20)         
           

