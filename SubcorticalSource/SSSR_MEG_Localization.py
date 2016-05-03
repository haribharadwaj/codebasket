import mne
import numpy as np
import os
import fnmatch
from anlffr import spectral
import pylab as pl
from scipy import io


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

# froot = '/home/hari/Documents/SubcorticalSource/'
froot = '/space/zazen/1/users/SubcorticalSource/'
saveRaw = True
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
        # raw.resample(raw.info['sfreq'] / 2.0)
        # Filter the data for SSRs
        raw.filter(l_freq=70., h_freq=None, filter_length='500ms', picks=None)

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
    params = dict(Fs=Fs, fpass=[100, 300], tapers=[1, 1], itc=0)

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
    figname_eeg = subj + '_' + condstem + '_eeg-results.pdf'
    pl.savefig(fpath + figname_eeg)
    y_eeg = y

    # Magnetometers
    pl.figure()
    for badname in epochs.info['bads']:
        bad_ind = epochs.info['ch_names'].index(badname)
        if bad_ind in mags:
            mags.remove(bad_ind)
    plv = np.zeros((len(mags), plv.shape[0]))
    for k, ch in enumerate(mags):
        y = x[:, [ch, ], :].transpose((1, 0, 2))
        print 'Computing PLV for channel', k
        plv[k, :], f = spectral.mtplv(y, params, verbose='DEBUG')
    f_AM = AMlist[cond - 1]
    ind_AM = np.argmin(np.abs(f - f_AM))
    best10Chans = plv[:, ind_AM].argsort()[:-10:-1]
    pl.plot(f, plv[best10Chans, :].mean(axis=0), linewidth=2)
    pl.xlabel('Frequency (Hz)', fontsize=16)
    pl.ylabel('Intertrial PLV', fontsize=16)
    pl.title('MEG Magnetometers', fontsize=16)
    pl.xlim([100, 300])
    pl.show()
    figname_mag = subj + '_' + condstem + '_mag-results.pdf'
    pl.savefig(fpath + figname_mag)
    lout = mne.find_layout(epochs.info)
    pos = lout.pos[mags]
    pl.figure()
    mne.viz.plot_topomap(plv[:, ind_AM], pos, sensors='ok', vmin=-0.01,
                         vmax=0.01)
    pl.show()
    figname_mag_topo = subj + '_' + condstem + '_mag_topo-results.pdf'
    pl.savefig(fpath + figname_mag_topo)
    y_mag = x[:, np.asarray(mags)[best10Chans], :].transpose((1, 0, 2))
    pos_mag = pos
    topo_mag = plv[:, ind_AM]

    doGrad = False
    if doGrad:
        # Gradiometers
        pl.figure()
        params = dict(Fs=Fs, fpass=[100, 300], tapers=[1, 1], itc=0)
        for badname in epochs.info['bads']:
            bad_ind = epochs.info['ch_names'].index(badname)
            if bad_ind in grads:
                grads.remove(bad_ind)

        plv = np.zeros((len(grads), plv.shape[1]))
        for k, ch in enumerate(grads):
            print 'Computing PLV for channel', k
            y = x[:, [ch, ], :].transpose((1, 0, 2))
            plv[k, :], f = spectral.mtplv(y, params, verbose='DEBUG')
        f_AM = AMlist[cond - 1]
        ind_AM = np.argmin(np.abs(f - f_AM))
        best10Chans = plv[:, ind_AM].argsort()[:-10:-1]
        pl.plot(f, plv[best10Chans, :].mean(axis=0), linewidth=2)
        pl.xlabel('Frequency (Hz)', fontsize=16)
        pl.ylabel('Intertrial PLV', fontsize=16)
        pl.title('MEG Gradiometers', fontsize=16)
        pl.xlim([100, 300])
        pl.show()
        figname_grad = subj + '_' + condstem + '_grad-results.pdf'
        pl.savefig(fpath + figname_grad)
        lout = mne.find_layout(epochs.info)
        pos = lout.pos[grads]
        pl.figure()
        mne.viz.plot_topomap(plv[:, ind_AM], pos, sensors='ok', vmin=-0.01,
                             vmax=0.01)
        pl.show()
        figname_grad_topo = (subj + '_' + condstem + '_grad_topo-results.pdf')
        pl.savefig(fpath + figname_grad_topo)
        y_grad = x[:, np.asarray(grads)[best10Chans], :].transpose((1, 0, 2))
        pos_grad = pos
        topo_grad = plv[:, ind_AM]
        saveDict = dict(Fs=Fs, params=params, y_mag=y_mag, y_grad=y_grad,
                        y_eeg=y_eeg, topo_mag=topo_mag, topo_grad=topo_grad,
                        f_AM=f_AM, t=times)
    else:
        saveDict = dict(Fs=Fs, params=params, y_mag=y_mag,
                        y_eeg=y_eeg, topo_mag=topo_mag,
                        f_AM=f_AM, t=times)
    save_res_name = subj + '_' + condstem + '_results.mat'
    io.savemat(fpath + save_res_name, saveDict)
    evoked = epochs.average()
    evoked.save(fpath + subj + '_' + condstem + '-ave.fif')
    cov = mne.cov.compute_covariance(epochs, tmin=None, tmax=0.002,
                                     verbose='DEBUG')
    cov.save(fpath + subj + '_' + condstem + '-cov.fif')
