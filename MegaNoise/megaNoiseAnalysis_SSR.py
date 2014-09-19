import mne
import numpy as np
from scipy import io
import os
import fnmatch
from anlffr import spectral
from anlffr.helper import biosemi2mne as bs
from anlffr.preproc import find_blinks
from mne.preprocessing.ssp import compute_proj_epochs
from mne.time_frequency.tfr import _induced_power as induced_power
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
froot = '/home/hari/Documents/PythonCodes/MegaNoise/'

lowLevel = False
if lowLevel:
    froot = froot + 'Level58dB/'

subjlist = ['I13', ]
nchans = 34
ch = [3, 4, 25, 26, 30, 31]  # Channels of interest
freqs = np.arange(5, 500, 2)  # define frequencies of interest
n_cycles = freqs / float(3)  # different number of cycle per frequency
n_cycles[freqs < 15] = 2

SSSR = False
ASSR25 = True  # Set false for ASSR43

for subj in subjlist:

    fpath = froot + subj + '/'

    # These are so that the generated files are organized better
    respath = fpath + 'RES/'

    if SSSR:
        condlist = [11, 12, 13, 14, 21, 22, 23, 24,
                    31, 32, 33, 34, 41, 42, 43, 44]
        condstem = 'allEvents'
    else:
        if ASSR25:
            condlist = [11, 13, 21, 23, 31, 33, 41, 43]
            condstem = 'ASSR25'
        else:
            condlist = [12, 14, 22, 24, 32, 34, 42, 44]
            condstem = 'ASSR43'

    print 'Running Subject', subj, 'Condition', condstem

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
            (raw, eves) = bs.importbdf(fpath + edfname, nchans=nchans,
                                       refchans=['EXG1', 'EXG2'],
                                       verbose='DEBUG')

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
            epochs = mne.Epochs(raw, eves, condlist, tmin=-0.1, proj=True,
                                tmax=1.0, baseline=(-0.1, 0.0),
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
    powtemp, plvtemp = induced_power(dat, sfreq=Fs, frequencies=freqs,
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
    bmin, bmax = -0.1, 0.0
    plvbline = plv[:, np.logical_and(times < bmax,
                                     times > bmin)].mean(axis=1)
    pl.figure()
    pl.imshow((plv.T - plvbline.T).T, vmin=0.02, vmax=0.2,
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

    #########################################################################
    # Fourier domain stuff
    pl.figure()
    params = dict(Fs=Fs, fpass=[5, 200], tapers=[1, 1], itc=1)
    y = x[:, :32, :].transpose((1, 0, 2))
    plv, f = spectral.mtcpca(y, params, verbose='DEBUG')
    pl.plot(f, plv**0.5, linewidth=2)
    pl.xlabel('Frequency (Hz)', fontsize=16)
    pl.ylabel('Intertrial PLV', fontsize=16)
    pl.show()
