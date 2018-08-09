from anlffr.helper import biosemi2mne as bs
import mne
import numpy as np
from anlffr import spectral
from scipy import io
import os
import fnmatch
import pylab as pl

# Adding Files and locations
froot = 'D:/DATA/EFR/'

subjlist = ['S050', ]

cond = 1
overwriteOld = False
for subj in subjlist:
    fpath = froot + subj + '/'

    # These are so that the generated files are organized better
    respath = fpath + 'RES/'

    print 'Running Subject', subj

    save_raw_name = subj + '_fmsweep_alltrial.mat'

    if os.path.isfile(respath + save_raw_name) and not overwriteOld:
        print 'Epoched data is already available on disk!'
        print 'Loading data from:', respath + save_raw_name
        x = io.loadmat(respath + save_raw_name)['x']
        fs = 4096.0
    else:
        bdfs = fnmatch.filter(os.listdir(fpath), subj + '*EFR_fmsweep*.bdf')
        print 'No pre-epoched data found, looking for BDF files'
        print 'Viola!', len(bdfs),  'files found!'

        for k, edfname in enumerate(bdfs):
            # Load data and read event channel
            (raw, eves) = bs.importbdf(fpath + edfname, nchans=36,
                                       refchans=['EXG1', 'EXG2'],
                                       mask=None)
            eves[:, 1] = np.mod(eves[:, 1], 256)
            eves[:, 2] = np.mod(eves[:, 2], 256)
            raw.set_channel_types({'EXG3': 'eeg', 'EXG4': 'eeg'})
            raw.info['bads'] += ['EXG3', 'A1', 'A2', 'A30', 'A7', 'A6',
                                 'A24', 'A28', 'A29', 'A3', 'A11', 'A15',
                                 'A16', 'A17', 'A10', 'A21', 'A20', 'A25',
                                 'EXG4']
            raw.drop_channels(raw.info['bads'])
            # Filter the data
            raw.filter(l_freq=130, h_freq=500)

            # raw.apply_proj()
            fs = raw.info['sfreq']

            # Epoching events of type
            epochs = mne.Epochs(
                raw, eves, cond, tmin=0., proj=False,
                tmax=3.995, baseline=(0., 3.995),
                reject=dict(eeg=100e-6))  # 200 regular, 50 strict

            xtemp = epochs.get_data()

            # Reshaping to the format needed by spectral.mtcpca() and
            # calling it
            if(xtemp.shape[0] > 0):
                xtemp = xtemp.transpose((1, 0, 2))
                xtemp = xtemp[:15, :, :]
                if(k == 0):
                    x = xtemp
                else:
                    x = np.concatenate((x, xtemp), axis=1)
            else:
                continue

    params = dict(Fs=fs, fpass=[5, 1000], tapers=[3, 5], Npairs=2000,
                  itc=0, nfft=fs*4)

    Ntrials = x.shape[1]

    print 'Running Mean Spectrum Estimation'
    (S, N, f) = spectral.mtspec(x, params, verbose=True)

    print 'Running channel by channel PLV Estimation'
    (plv, f) = spectral.mtplv(x, params, verbose=True)

    # Saving Results
    res = dict(plv=plv, f=f, S=S, N=N, Ntrials=Ntrials)

    save_name = subj + '_fmsweep_results.mat'

    if (not os.path.isdir(respath)):
        os.mkdir(respath)
    io.savemat(respath + save_name, res)

    if not os.path.isfile(respath + save_raw_name):
        io.savemat(respath + save_raw_name, dict(x=x, subj=subj))

    plotStuff = True
    if plotStuff:
        pl.figure()
        pl.plot(f, S[-1, ])
        pl.hold('on')
        pl.plot(f, N[-1, ])
        pl.xlabel('Frequency (Hz)', fontsize=16)
        pl.ylabel('Power', fontsize=16)
        ax = pl.gca()
        ax.set_xlim([100, 500])
        ax.tick_params(labelsize=16)
        pl.show()
