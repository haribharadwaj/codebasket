import mne
import numpy as np
import os
import fnmatch
from anlffr import spectral
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
froot = '/autofs/cluster/transcend/hari/ASSRnew/'
saveResults = True
subjlist = ['manny', ]
ch = range(1, 307)  # Channels of interest
mags = range(2, 306, 3)
grads = range(0, 306, 3) + range(1, 306, 3)
freqs = np.arange(5, 500, 2)  # define frequencies of interest
n_cycles = freqs / float(3)  # different number of cycle per frequency
n_cycles[freqs < 15] = 2

SSSR = False
ASSR25 = False  # Set false for ASSR43
sss = False
eeg = False
for subj in subjlist:

    fpath = froot + subj + '/'
    if sss:
        ssstag = '_sss'
    else:
        ssstag = ''

    if SSSR:
        condlist = range(1, 13)
        condstem = 'allEvents'
        l_freq = 70
    else:
        l_freq = 2.0
        if ASSR25:
            condlist = [1, 3, 5, 7, 9, 11]
            condstem = 'ASSR25'
        else:
            condlist = [2, 4, 6, 8, 10, 12]
            condstem = 'ASSR43'

    print 'Running Subject', subj, 'Condition', condstem

    save_raw_name = subj + '_' + condstem + ssstag + '-epo.fif'

    if os.path.isfile(fpath + save_raw_name):
        preEpoched = True
        print 'Epoched data is already available on disk!'
        print 'Loading data from:', fpath + save_raw_name
        epochs = mne.read_epochs(fpath + save_raw_name, verbose='DEBUG')
        Fs = epochs.info['sfreq']
        x = epochs.get_data()
        times = epochs.times
    else:
        preEpoched = False
        fifs = fnmatch.filter(os.listdir(fpath), subj + '*raw' +
                              ssstag + '.fif')
        print 'No pre-epoched data found, looking for raw files'
        print 'Viola!', len(fifs),  'files found!'
        for k, fif in enumerate(fifs):
            fifs[k] = fpath + fif
        # Load data and read event channel
        raw = mne.io.Raw(fifs, preload=True, add_eeg_ref=False)
        eves = mne.find_events(raw, stim_channel='STI101',
                               shortest_event=1)
        if not sss:
            raw.info['bads'] += ['MEG2033', 'MEG0442', 'MEG2343', 'MEG1643',
                                 'MEG1211', 'MEG2522', 'MEG0731']
        if eeg:
            raw.info['bads'] += ['EEG004', 'EEG038', 'EEG040', 'EEG067']
        # Filter the data for SSRs
        raw.filter(l_freq=l_freq, h_freq=144, l_trans_bandwidth=0.15,
                   picks=np.arange(0, 306, 1))

        if not SSSR:
            # raw.resample(sfreq=1000.0, n_jobs=4, verbose='DEBUG')
            # SSP for blinks
            blinks = find_blinks(raw, ch_name='EOG062')
            epochs_blinks = mne.Epochs(raw, blinks, 998, tmin=-0.25,
                                       tmax=0.25, proj=True,
                                       baseline=(-0.25, 0),
                                       reject=dict(grad=8000e-13,
                                                   mag=8e-12))
            blink_projs = compute_proj_epochs(epochs_blinks, n_grad=2,
                                              n_mag=2, n_eeg=0,
                                              verbose='DEBUG')
            raw.add_proj(blink_projs)

            # SSP for cardiac
            qrs = find_blinks(raw, ch_name='ECG063', h_freq=100.0,
                              event_id=999)
            epochs_qrs = mne.Epochs(raw, qrs, 999, tmin=-0.1,
                                    tmax=0.1, proj=True,
                                    baseline=(-0.1, 0),
                                    reject=dict(grad=8000e-13,
                                                mag=8e-12))
            qrs_projs = compute_proj_epochs(epochs_qrs, n_grad=2,
                                            n_mag=2, n_eeg=0,
                                            verbose='DEBUG')
            raw.add_proj(qrs_projs)
            useProj = True
        else:
            useProj = False

        # Epoching events of type
        epochs = mne.Epochs(raw, eves, condlist, tmin=-0.1, proj=useProj,
                            tmax=1.3, baseline=(-0.1, 0.),
                            reject=dict(grad=5000e-13, mag=4e-12))

        x = epochs.get_data()
        if saveResults:
            epochs.save(fpath + save_raw_name)

    # Calculate power, plv
    if not preEpoched:
        Fs = raw.info['sfreq']
        times = epochs.times

    calculateTF = False
    if calculateTF:
        plv = np.zeros((len(freqs), len(times)))
        tfspec = np.zeros((len(freqs), len(times)))
        dat = x[:, ch, :].mean(axis=1, keepdims=True)
        powtemp, plvtemp = induced_power(dat, sfreq=Fs, frequencies=freqs,
                                         n_cycles=n_cycles, zero_mean=True)
        plv = plvtemp.squeeze()
        tfspec = pow2db(powtemp.squeeze())

        ######################################################################
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
    params = dict(Fs=Fs, fpass=[5, 140], tapers=[1, 1], itc=0)
    for badname in epochs.info['bads']:
        bad_ind = epochs.info['ch_names'].index(badname)
        if bad_ind in mags:
            mags.remove(bad_ind)

    y = x[:, mags, :].transpose((1, 0, 2))
    plv, f = spectral.mtplv(y, params, verbose='DEBUG')
    pl.plot(f, plv.T, linewidth=2)
    pl.xlabel('Frequency (Hz)', fontsize=16)
    pl.ylabel('Intertrial PLV', fontsize=16)
    pl.title('MEG Magnetometers', fontsize=16)
    if SSSR:
        pl.xlim([70, 140])
    else:
        pl.xlim([5, 140])
    pl.show()
    lout = mne.find_layout(epochs.info)
    pos = lout.pos[mags]
    if SSSR:
        f_AM = 107.0
    else:
        if ASSR25:
            f_AM = 25.0
        else:
            f_AM = 43.0
    ind_AM = np.argmin(np.abs(f - f_AM))
    pl.figure()
    mne.viz.plot_topomap(plv[:, ind_AM], pos, sensors='ok', vmin=-0.07,
                         vmax=0.07, show_names=True)
    pl.show()

    if eeg:
        eeg = []
        for k, ch in enumerate(raw.info['ch_names']):
            if 'EEG' in ch:
                eeg += [k, ]
        yy = x[:, eeg, :].transpose((1, 0, 2))
        pl.figure()
        plv, f = spectral.mtplv(yy, params, verbose='DEBUG')
        pl.plot(f, plv.T, linewidth=2)
        pl.xlabel('Frequency (Hz)', fontsize=16)
        pl.ylabel('Intertrial PLV', fontsize=16)
        pl.title('EEG', fontsize=16)
        if SSSR:
            pl.xlim([70, 140])
        else:
            pl.xlim([5, 140])
        pl.show()
        louteeg = mne.layouts.make_eeg_layout(raw.info,
                                              exclude=['EEG073', 'EEG074',
                                                       'EEG097', 'EEG098'])
        pl.figure()
        mne.viz.plot_topomap(plv[:-4, ind_AM], louteeg.pos, sensors='ok',
                             vmin=-0.2, vmax=0.2)
        pl.show()
