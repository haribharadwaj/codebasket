from anlffr.helper import biosemi2mne as bs
import mne
import numpy as np
from anlffr import spectral
from scipy import io
import os
import fnmatch
import pylab as pl

# Adding Files and locations
# froot = '/home/hari/Documents/PythonCodes/MEMR/'
froot = '/autofs/cluster/transcend/hari/MEMR/'


subjlist = ['I41', ]
ear = 'right'

for subj in subjlist:

    fpath = froot + subj + '/'

    # These are so that the generated files are organized better
    respath = fpath + 'RES/'

    cond = 1
    condstem = '_223Hz_EFR'
    print 'Running Subject', subj, 'Condition', cond

    save_raw_name = subj + '_' + ear + condstem + '_alltrial.mat'

    if os.path.isfile(respath + save_raw_name):
        print 'Epoched data is already available on disk!'
        print 'Loading data from:', respath + save_raw_name
        x = io.loadmat(respath + save_raw_name)['x']
        fs = 4096.0
    else:
        bdfs = fnmatch.filter(os.listdir(fpath), subj +
                              '*' + ear + '*EFR*.bdf')
        if len(bdfs) == 0:
            bdfs = fnmatch.filter(os.listdir(fpath), subj +
                                  '*EFR*' + ear + '*.bdf')

        print 'No pre-epoched data found, looking for BDF files'
        print 'Viola!', len(bdfs),  'files found!'

        for k, edfname in enumerate(bdfs):
            # Load data and read event channel
            (raw, eves) = bs.importbdf(fpath + edfname, nchans=35,
                                       refchans=['EXG1', 'EXG2'])

            raw.info['bads'] += ['EXG3', 'A6', 'A7', 'A24', 'A25']
            # Filter the data
            raw.filter(
                l_freq=70, h_freq=1500, picks=np.arange(0, 35, 1))

            # raw.apply_proj()
            fs = raw.info['sfreq']

            # Epoching events of type
            epochs = mne.Epochs(
                raw, eves, cond, tmin=-0.025, proj=False,
                tmax=1.025, baseline=(-0.025, 0.0),
                reject=dict(eeg=200e-6))

            xtemp = epochs.get_data()

            # Reshaping to the format needed by spectral.mtcpca() and
            # calling it
            if(xtemp.shape[0] > 0):
                xtemp = xtemp.transpose((1, 0, 2))
                xtemp = xtemp[0:32, :, :]
                if(k == 0):
                    x = xtemp
                else:
                    x = np.concatenate((x, xtemp), axis=1)
            else:
                continue

    nPerDraw = 400
    nDraws = 100
    params = dict(Fs=fs, fpass=[5, 1000], tapers=[1, 1], Npairs=2000,
                  itc=1, nfft=8192)

    Ntrials = x.shape[1]

    print 'Running Mean Spectrum Estimation'
    (S, N, f) = spectral.mtspec(x, params, verbose=True)

    print 'Running CPCA PLV Estimation'
    (cplv, f) = spectral.mtcpca(x, params, verbose=True)

    print 'Running channel by channel PLV Estimation'
    (plv, f) = spectral.mtplv(x, params, verbose=True)

    print 'Running CPCA Power Estimation'
    (cpow, f) = spectral.mtcspec(x, params, verbose=True)

    print 'Running raw spectrum estimation'
    (Sraw, f) = spectral.mtspecraw(x, params, verbose=True)

    # Saving Results
    res = dict(cpow=cpow, plv=plv, cplv=cplv, Sraw=Sraw,
               f=f, S=S, N=N, Ntrials=Ntrials)

    save_name = subj + '_' + ear + condstem + '_results.mat'

    if (not os.path.isdir(respath)):
        os.mkdir(respath)
    io.savemat(respath + save_name, res)

    if not os.path.isfile(respath + save_raw_name):
        io.savemat(respath + save_raw_name, dict(x=x, subj=subj))

    pl.figure()
    pl.plot(f, cplv, linewidth=2)
    pl.xlabel('Frequency (Hz)', fontsize=16)
    pl.ylabel('Phase Locking', fontsize=16)
    ax = pl.gca()
    ax.set_xlim([200, 500])
    ax.tick_params(labelsize=16)
    pl.show()
