from anlffr.helper import biosemi2mne as bs
import mne
import numpy as np
import os
import fnmatch
from scipy.signal import savgol_filter as sg
from scipy.io import savemat
from anlffr.spectral import mtplv
import pylab as pl


# Setup bayesian-weighted averaging
def bayesave(x, trialdim=0, timedim=1, method='mean', smoothtrials=19):
    ntrials = x.shape[trialdim]
    if method == 'mean':
        summary = x.mean(axis=trialdim) * (ntrials) ** 0.5
    else:
        summary = np.median(x, axis=trialdim) * (ntrials) ** 0.5

    ntime = x.shape[timedim]
    wts = 1 / np.var(x - summary, axis=timedim, keepdims=True)
    wts = sg(wts, smoothtrials, 3, axis=trialdim)  # Smooth the variance
    normFactor = wts.sum()
    wts = wts.repeat(ntime, axis=timedim)
    ave = (x * wts).sum(axis=trialdim) / normFactor
    return ave


# NOTE: THIS IS EARPLUG DATA CODE.. NOT REGULAR CENTRAL GAIN CODE!!!!
# Adding Files and locations
froot = '/Users/HMB105/Library/CloudStorage/Dropbox/Data_Testing/Earplug_ABRs/'

# NOTE: THIS IS EARPLUG DATA CODE.. NOT REGULAR CENTRAL GAIN CODE!!!!
# Done: ['E002', 'E004', 'E005', 'E006', 'E012', 'E014','E022', 'E025', 'Hari'
subjlist = ['E007']

earlist = ['L', 'R']

visitlist = [1, 2]
for visit in visitlist:
    for subj in subjlist:
        for ear in earlist:
            if ear == 'L':
                conds = [[6, 10], [9, 12]]
                names = ['_L_moderate', '_L_loud']
            else:
                conds = [[96, 160], [144, 192]]
                names = ['_R_moderate', '_R_loud']
    
            print(f'Running Subject {subj}, {ear} ear')
            for ind, cond in enumerate(conds):
                name = names[ind]
                print(f'Doing condition {cond}')
                fpath = froot + '/' + subj + '/'
                
                bdfs = fnmatch.filter(os.listdir(fpath), subj + f'_ABR_{visit}*.bdf')
    
                if len(bdfs) >= 1:
                    for k, bdf in enumerate(bdfs):
                        edfname = fpath + bdf
                        # Load data and read event channel
                        extrachans = [u'GSR1', u'GSR2', u'Erg1', u'Erg2', u'Resp',
                                      u'Plet', u'Temp']
                        raw, eves = bs.importbdf(edfname, nchans=36,
                                                 extrachans=extrachans)
                        raw.set_channel_types({'EXG3': 'eeg', 'EXG4': 'eeg'})
    
                        # Pick channels to not include in epoch rejection
                        raw.info['bads'] += ['EXG3', 'EXG4', 'A1', 'A30',
                                             'A15', 'A16']
                        # Filter the data
                        raw.filter(l_freq=30., h_freq=3000, picks=np.arange(36))
    
                        # Epoch the data
                        tmin, tmax = -0.002, 0.08
                        bmin, bmax = -0.002, 0.0016
                        rejthresh = 100e-6
    
                        # Check if condition has events in the file:
                        if ((np.equal(eves[:, 2], cond[0]).sum() == 0) or
                            (np.equal(eves[:, 2], cond[0]).sum() == 0)):
                            continue
    
                        epochs = mne.Epochs(raw, eves, cond, tmin=tmin, proj=False,
                                            tmax=tmax, baseline=(bmin, bmax),
                                            picks=np.arange(36),
                                            reject=dict(eeg=rejthresh),
                                            verbose='WARNING')
                        xtemp = epochs.get_data()
                        t = epochs.times * 1e3 - 1.6  # Adjust for delay and use ms
                        # Reshaping so that channels is first
                        if(xtemp.shape[0] > 0):
                            xtemp = xtemp.transpose((1, 0, 2))
                            if ('x' not in locals()) or (k == 0):
                                x = xtemp
                            else:
                                x = np.concatenate((x, xtemp), axis=1)
                        else:
                            continue
                else:
                    RuntimeError('No BDF files found!!')
    
                # Average data
                # goods = [1, 5, 6, 27, 28, 30, 31]
                
                # Pitt channels
                goods = [28, 3, 30, 26, 4, 25, 7, 31, 22, 9, 8, 21, 11, 12, 18]
                
                if ear == 'L':
                    refchan = 34
                else:
                    refchan = 35
    
                y = x[goods, :, :].mean(axis=0) - x[refchan, :, :]
    
                params = dict(Fs=raw.info['sfreq'], tapers=[4, 7],
                              fpass=[1, 4000], itc=0)
                plv, f = mtplv(y, params)
                plv_mean = sg(plv, 5, 3)
    
                tdresp = bayesave(x[goods, :, :].mean(axis=0)) * 1e6  # Microvolts

    
                # Make dictionary and save
                mdict = dict(t=t, x=tdresp, f=f, plv=plv_mean)
                savepath = froot + '/CentralGainResults/'
                savename = subj + name + f'_{visit}_CentralGain.mat'
                savemat(savepath + savename, mdict)
    
                pl.plot(f, plv, linewidth=2)
                pl.xlabel('Modulation Frequency (Hz)', fontsize=16)
                pl.xticks(fontsize=14)
                pl.yticks(fontsize=14)
                pl.xscale('log')
                pl.xticks([16, 32, 64, 128, 256, 512, 1024, 2048],
                          labels=['16', '32', '64', '128', '256', '512', '1024', '2048'],
                          fontsize=14)
                pl.xlim((16, 2048))
                pl.grid()


