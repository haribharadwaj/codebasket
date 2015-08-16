import mne
import numpy as np
import os
import fnmatch
from anlffr import spectral
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
froot = '/autofs/cluster/transcend/hari/MEMR/'
# froot = '/home/hari/Documents/SubcorticalSource/'
saveRaw = False
subjlist = ['HB', ]

ch = range(1, 307)  # Channels of interest
mags = range(2, 306, 3)
grads = range(0, 306, 3) + range(1, 306, 3)
eeg = [306, ]
freqs = np.arange(5, 500, 2)  # define frequencies of interest
n_cycles = freqs / float(3)  # different number of cycle per frequency
n_cycles[freqs < 15] = 2

cond = 3
AMlist = [163, 193, 223, 253, 283]
for subj in subjlist:
    fpath = froot + subj + '/'
    condstem = str(AMlist[cond - 1]) + 'Hz'

    print 'Running Subject', subj, 'Condition', condstem

    save_raw_name = subj + '_' + condstem + '-epo.fif'

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
        fifs = fnmatch.filter(os.listdir(fpath), subj + '_ASSR_*' +
                              str(cond) + '_raw.fif')
        print 'No pre-epoched data found, looking for raw files'
        print 'Viola!', len(fifs),  'files found!'
        for k, fif in enumerate(fifs):
            fifs[k] = fpath + fif
        # Load data and read event channel
        raw = mne.io.Raw(fifs, preload=True, add_eeg_ref=False)
        eves = mne.find_events(raw, stim_channel='STI101',
                               shortest_event=1)

        raw.set_channel_types({'EMG061': 'eeg'})
        raw.info['bads'] += ['MEG0223', 'MEG1623']
        raw.resample(raw.info['sfreq'] / 2.0)
        # Filter the data for SSRs
        # raw.filter(l_freq=70., h_freq=300., picks=None)

        # Epoching events of type
        epochs = mne.Epochs(raw, eves, cond, tmin=-0.05, proj=False,
                            tmax=1.25, baseline=(-0.05, 0.),
                            reject=dict(grad=5000e-13, mag=4e-12))

        x = epochs.get_data()
        if saveRaw:
            epochs.save(fpath + save_raw_name)

    # Calculate power, plv
    if not preEpoched:
        Fs = raw.info['sfreq']
        times = epochs.times

    #########################################################################
    # Fourier domain stuff
    pl.figure()
    params = dict(Fs=Fs, fpass=[100, 300], tapers=[1, 1], itc=0)
    for badname in epochs.info['bads']:
        bad_ind = epochs.info['ch_names'].index(badname)
        if bad_ind in mags:
            mags.remove(bad_ind)

    y = x[:, mags, :].transpose((1, 0, 2))
    plv, f = spectral.mtplv(y, params, verbose='DEBUG')
    pl.plot(f, plv.mean(axis=0), linewidth=2)
    pl.xlabel('Frequency (Hz)', fontsize=16)
    pl.ylabel('Intertrial PLV', fontsize=16)
    pl.title('MEG Magnetometers', fontsize=16)
    pl.xlim([100, 300])
    pl.show()
    figname_mag = subj + '_' + condstem + '_mag-results.pdf'
    pl.savefig(fpath + figname_mag)
    lout = mne.find_layout(epochs.info)
    pos = lout.pos[mags]
    f_AM = AMlist[cond - 1]
    ind_AM = np.argmin(np.abs(f - f_AM))
    pl.figure()
    mne.viz.plot_topomap(plv[:, ind_AM], pos, sensors='ok', vmin=-0.01,
                         vmax=0.01)
    pl.show()
    figname_mag_topo = subj + '_' + condstem + '_mag_topo-results.pdf'
    pl.savefig(fpath + figname_mag_topo)

    # Gradiometers
    pl.figure()
    params = dict(Fs=Fs, fpass=[100, 300], tapers=[1, 1], itc=0)
    for badname in epochs.info['bads']:
        bad_ind = epochs.info['ch_names'].index(badname)
        if bad_ind in grads:
            grads.remove(bad_ind)

    y = x[:, grads, :].transpose((1, 0, 2))
    plv, f = spectral.mtplv(y, params, verbose='DEBUG')
    pl.plot(f, plv.mean(axis=0), linewidth=2)
    pl.xlabel('Frequency (Hz)', fontsize=16)
    pl.ylabel('Intertrial PLV', fontsize=16)
    pl.title('MEG Gradiometers', fontsize=16)
    pl.xlim([100, 300])
    pl.show()
    figname_grad = subj + '_' + condstem + '_grad-results.pdf'
    pl.savefig(fpath + figname_grad)
    lout = mne.find_layout(epochs.info)
    pos = lout.pos[grads]
    f_AM = AMlist[cond - 1]
    ind_AM = np.argmin(np.abs(f - f_AM))
    pl.figure()
    mne.viz.plot_topomap(plv[:, ind_AM], pos, sensors='ok', vmin=-0.01,
                         vmax=0.01)
    pl.show()
    figname_grad_topo = (subj + '_' + condstem + '_grad_topo-results.pdf')
    pl.savefig(fpath + figname_grad_topo)

    # EEG
    pl.figure()
    params = dict(Fs=Fs, fpass=[100, 300], tapers=[1, 1], itc=0)
    y = x[:, eeg, :].transpose((1, 0, 2))
    plv, f = spectral.mtplv(y, params, verbose='DEBUG')
    pl.plot(f, plv.T, linewidth=2)
    pl.xlabel('Frequency (Hz)', fontsize=16)
    pl.ylabel('Intertrial PLV', fontsize=16)
    pl.title('EEG', fontsize=16)
    pl.xlim([100, 300])
    pl.show()
    figname_grad = subj + '_' + condstem + '_eeg-results.pdf'
    pl.savefig(fpath + figname_grad)
