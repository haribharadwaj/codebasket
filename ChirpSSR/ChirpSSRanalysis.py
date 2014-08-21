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


def pow2db(x):
    """ Converts *power* to decibels

    Parameters
    ----------

    x - Input in linear units

    Returns
    -------

    y - Equivalend in decibel units
    """

    y = 10*np.log10(x)
    return y

# Adding Files and locations
froot = '/home/hari/Documents/PythonCodes/ChirpSSR/'

lowLevel = True
if lowLevel:
    froot = froot + 'Level58dB/'

subjlist = ['I41', ]
ch = [3, 4, 25, 26, 30, 31]  # Channels of interest
freqs = np.arange(2, 500, 2)  # define frequencies of interest
n_cycles = freqs / float(5)  # different number of cycle per frequency
n_cycles[freqs < 15] = 2

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
            preEpoched = True
            print 'Epoched data is already available on disk!'
            print 'Loading data from:', respath + save_raw_name
            dat = io.loadmat(respath + save_raw_name)
            Fs = dat['Fs'][0, 0]
            x = dat['x']
            times = dat['times'].squeeze()
        else:
            preEpoched = False
            bdfs = fnmatch.filter(os.listdir(fpath), subj + '*.bdf')
            print 'No pre-epoched data found, looking for BDF files'
            print 'Viola!', len(bdfs),  'files found!'

            for k, edfname in enumerate(bdfs):
                # Load data and read event channel
                (raw, eves) = bs.importbdf(fpath + edfname, nchans=35,
                                           refchans=['EXG1', 'EXG2'])

                raw.info['bads'] += ['EXG3', ]
                # Filter the data for ERPs
                raw.filter(l_freq=0.5, h_freq=500, l_trans_bandwidth=0.15,
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
                raw.del_proj(0)  # Removing average reference projection
                raw.add_proj(blink_projs)

                # Epoching events of type
                epochs = mne.Epochs(raw, eves, cond, tmin=-0.1, proj=True,
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
    if not preEpoched:
        Fs = raw.info['sfreq']
        times = epochs.times
    plv = np.zeros((len(freqs), len(times)))
    tfspec = np.zeros((len(freqs), len(times)))
    dat = x[:, ch, :].mean(axis=1, keepdims=True)
    powtemp, plvtemp = induced_power(dat, Fs=Fs, frequencies=freqs,
                                     n_cycles=n_cycles, zero_mean=True)
    plv = plvtemp.squeeze()
    tfspec = pow2db(powtemp.squeeze())

    # Save results to the RES directory
    savedict = dict(tfspec=tfspec, plv=plv, times=times, freqs=freqs)
    save_name = subj + '_plv_inducedpow_2Hz_500Hz.mat'
    print 'Saving data for subject', subj
    if (not os.path.isdir(respath)):
        os.mkdir(respath)
    io.savemat(respath + save_name, savedict)

    if not os.path.isfile(respath + save_raw_name):
        io.savemat(respath + save_raw_name, dict(x=x, subj=subj, Fs=Fs,
                                                 times=times))

    ##########################################################################
    # View time-frequency plots
    # PLV
    pl.close('all')
    t0 = 0.0
    pl.figure()
    pl.imshow(plv, vmin=0.05, vmax=0.2,
              extent=[times[0] - t0, times[-1] - t0, freqs[0], freqs[-1]],
              aspect='auto', origin='lower')
    pl.xlabel('Time (s)')
    pl.ylabel('Frequency (Hz)')
    pl.title('Phase Locking')
    pl.colorbar()
    pl.show()

    # Induced power
    pl.figure()
    bmin, bmax = -0.1, 0.0
    powbline = tfspec[:, np.logical_and(times < bmax,
                                        times > bmin)].mean(axis=1)
    pl.imshow((tfspec.T-powbline.T).T, vmin=-3.0, vmax=3.0,
              extent=[times[0] - t0, times[-1] - t0, freqs[0], freqs[-1]],
              aspect='auto', origin='lower')
    pl.xlabel('Time (s)')
    pl.ylabel('Frequency (Hz)')
    pl.title('Power (dB re: baseline)')
    pl.colorbar()
    pl.show()
