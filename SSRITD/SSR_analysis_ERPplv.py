import mne
import numpy as np
# from anlffr import spectral
from scipy import io
import os
import fnmatch
# import pylab as pl

from anlffr.helper import biosemi2mne as bs
from anlffr.preproc import find_blinks
# from scipy.signal import detrend
from mne.preprocessing.ssp import compute_proj_epochs
from mne.time_frequency import induced_power

# Adding Files and locations
froot = '/home/hari/Documents/MATLAB/SSRITD/EEG/'

# List of files stems, each will be appended by run number
# Could be different for edf, fiff, eve etc.
# Use list [] and enumerate over if filenames are weird

subjlist = [
    'I02', 'I03', 'I05', 'I06', 'I08', 'I09', 'I11', 'I13', 'I14', 'I15',
    'I17', 'I18', 'I19', 'I20', 'I25', 'I26', 'I29', 'I30', 'I33', 'I35',
    'I36', 'I37']
#subjlist = ['I13', ]

ch = [3, 4, 25, 26, 30, 31]  # Channels of interest
plv_allsubj = []
pow_allsubj = []
freqs = np.arange(5, 50, 2)  # define frequencies of interest
n_cycles = freqs / float(5)  # different number of cycle per frequency

for k, subj in enumerate(subjlist):
    fpath = froot + subj + '/'

    # These are so that the generated files are organized better
    respath = fpath + 'RES/'

    print 'Running Subject', subj

    bdfs = fnmatch.filter(os.listdir(fpath), subj + '*.bdf')

    if len(bdfs) > 1:
        print '***WARNING!! Multiple files found! Cannot continue!***'
    else:
        print 'Viola! Data files found!'
        edfname = bdfs[0]

    # Load data and read event channel
    (raw, eves) = bs.importbdf(fpath + edfname)

    # Filter the data for ERPs
    raw.filter(l_freq=0.5, h_freq=200, l_trans_bandwidth=0.15,
               picks=np.arange(0, 32, 1))
    #raw.apply_function(detrend, picks=np.arange(0, 32, 1), dtype=None,
    #                  n_jobs=1, verbose=True)

    # SSP for blinks
    blinks = find_blinks(raw)
    epochs_blinks = mne.Epochs(raw, blinks, 998, tmin=-0.25, tmax=0.25,
                               proj=True, baseline=(-0.25, 0),
                               reject=dict(eeg=500e-6))
    blink_projs = compute_proj_epochs(epochs_blinks, n_grad=0, n_mag=0,
                                      n_eeg=2, verbose='DEBUG')
    raw.add_proj(blink_projs)

    # Epoching and averaging
    conds = [8, 9, 10, 11, 12]  # 50, 100, 200, 400, 800 microseconds
    condstr = map(str, conds)

    epochs = mne.Epochs(raw, eves, conds, tmin=-0.5, proj=True,
                        tmax=2.5, baseline=(-0.5, 0.0),
                        reject=dict(eeg=200e-6))

    # Calculate power, plv
    Fs = raw.info['sfreq']
    times = epochs.times
    plv = np.zeros((len(conds), len(freqs), len(times)))
    tfspec = np.zeros((len(conds), len(freqs), len(times)))

    condorder = np.zeros(len(conds))  # To save the order in which they come
    for condind, cond in enumerate(condstr):
        dat = epochs[cond].get_data()[:, ch, :].mean(axis=1, keepdims=True)
        powtemp, plvtemp = induced_power(dat, Fs=Fs, frequencies=freqs,
                                         n_cycles=n_cycles, zero_mean=True)
        plv[condind, :, :] = plvtemp.squeeze()
        tfspec[condind, :, :] = powtemp.squeeze()

    # Save results to the RES directory
    savedict = dict(tfspec=tfspec, plv=plv, times=times, conds=conds,
                    freqs=freqs)
    save_name = subj + '_plv_inducedpow.mat'
    print 'Saving data for subject', subj
    if (not os.path.isdir(respath)):
        os.mkdir(respath)
    io.savemat(respath + save_name, savedict)
