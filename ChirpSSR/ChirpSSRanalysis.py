import mne
import numpy as np
from scipy import io
import os
import fnmatch
from anlffr.helper import biosemi2mne as bs
from anlffr.preproc import find_blinks
from mne.preprocessing.ssp import compute_proj_epochs
from mne.time_frequency import induced_power
import pylab as pl


# Adding Files and locations
froot = '/home/hari/Documents/PythonCodes/ChirpSSR/'

subjlist = ['I41', ]
ch = [3, 4, 25, 26, 30, 31]  # Channels of interest
freqs = np.arange(2, 500, 2)  # define frequencies of interest
n_cycles = freqs / float(5)  # different number of cycle per frequency
n_cycles[freqs < 5] = 1

for subj in subjlist:

    fpath = froot + subj + '/'

    # These are so that the generated files are organized better
    respath = fpath + 'RES/'

    condlist = np.arange(1, 2)
    condstemlist = ['chirp20to400Hz', ]
    for condind, cond in enumerate(condlist):
        condstem = condstemlist[condind]
        print 'Running Subject', subj, 'Condition', condind

        save_raw_name = subj + '_' + condstem + '_alltrial.mat'

        if os.path.isfile(respath + save_raw_name):
            print 'Epoched data is already available on disk!'
            print 'Loading data from:', respath + save_raw_name
            x = io.loadmat(respath + save_raw_name)['x']
        else:
            bdfs = fnmatch.filter(os.listdir(fpath), subj + '*.bdf')
            print 'No pre-epoched data found, looking for BDF files'
            print 'Viola!', len(bdfs),  'files found!'

            for k, edfname in enumerate(bdfs):
                # Load data and read event channel
                (raw, eves) = bs.importbdf(fpath + edfname, nchans=35,
                                           refchans=['EXG1', 'EXG2'])

                raw.info['bads'] += ['EXG3', ]
                # Filter the data for ERPs
                raw.filter(l_freq=0.5, h_freq=200, l_trans_bandwidth=0.15,
                           picks=np.arange(0, 32, 1))

                # raw.apply_proj()
                fs = raw.info['sfreq']
                # SSP for blinks
                blinks = find_blinks(raw)
                epochs_blinks = mne.Epochs(raw, blinks, 998, tmin=-0.25,
                                           tmax=0.25, proj=True,
                                           baseline=(-0.25, 0),
                                           reject=dict(eeg=500e-6))
                blink_projs = compute_proj_epochs(epochs_blinks, n_grad=0,
                                                  n_mag=0, n_eeg=2,
                                                  verbose='DEBUG')
                raw.add_proj(blink_projs)

                # Epoching events of type
                epochs = mne.Epochs(
                    raw, eves, cond, tmin=-0.1, proj=False,
                    tmax=2.1, baseline=(-0.1, 0.0),
                    reject = dict(eeg=150e-6))

                xtemp = epochs.get_data()

                # Reshaping to the format needed by spectral.mtcpca() and
                # calling it
                if(xtemp.shape[0] > 0):
                    if(k == 0):
                        x = xtemp
                    else:
                        x = np.concatenate((x, xtemp), axis=0)
                else:
                    continue

    # Calculate power, plv
    Fs = raw.info['sfreq']
    times = epochs.times
    plv = np.zeros((len(freqs), len(times)))
    tfspec = np.zeros((len(freqs), len(times)))
    dat = x[:, ch, :].mean(axis=1, keepdims=True)
    powtemp, plvtemp = induced_power(dat, Fs=Fs, frequencies=freqs,
                                     n_cycles=n_cycles, zero_mean=True)
    plv = plvtemp.squeeze()
    tfspec = powtemp.squeeze()

    # Save results to the RES directory
    savedict = dict(tfspec=tfspec, plv=plv, times=times, freqs=freqs)
    save_name = subj + '_plv_inducedpow_2Hz_500Hz.mat'
    print 'Saving data for subject', subj
    if (not os.path.isdir(respath)):
        os.mkdir(respath)
    io.savemat(respath + save_name, savedict)

    if not os.path.isfile(respath + save_raw_name):
        io.savemat(respath + save_raw_name, dict(x=x, subj=subj, Fs=Fs))

    ##########################################################################
    # View time-frequency plots

    pl.close('all')
    t0 = 0.0
    pl.figure()
    pl.imshow(plv, vmin=0.1, vmax=0.20,
              extent=[times[0] - t0, times[-1] - t0, freqs[0], freqs[-1]],
              aspect='auto', origin='lower')
    pl.xlabel('Time (s)')
    pl.ylabel('Frequency (Hz)')
    pl.title('Phase Locking')
    pl.colorbar()
    pl.show()
