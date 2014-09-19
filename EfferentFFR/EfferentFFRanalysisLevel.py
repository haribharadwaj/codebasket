from anlffr.helper import biosemi2mne as bs
import mne
import numpy as np
from anlffr import spectral
from scipy import io
import os
import fnmatch

# Adding Files and locations
froot = '/home/hari/Documents/PythonCodes/EfferentFFR/LevelControl/'

carr100 = True
noise = True

if carr100:
    if noise:
        froot = froot + 'Noisecarrier100/'
    else:
        froot = froot + 'PTCarr100/'

subjlist = ['I14', ]

for subj in subjlist:

    fpath = froot + subj + '/'

    # These are so that the generated files are organized better
    respath = fpath + 'RES/'

    condlist = np.arange(1, 4)
    condstemlist = ['70dBSPL', '64dBSPL', '58dBSPL']
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

                raw.info['bads'] += ['EXG3', 'A6', 'A7', 'A24', 'A21']
                # Filter the data
                raw.filter(
                    l_freq=70, h_freq=1500, picks=np.arange(0, 35, 1))

                # raw.apply_proj()
                fs = raw.info['sfreq']

                # Epoching events of type
                epochs = mne.Epochs(
                    raw, eves, cond, tmin=-0.025, proj=False,
                    tmax=0.325, baseline=(-0.025, 0.0),
                    reject = dict(eeg=150e-6))

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
        fs = 4096.0
        params = dict(Fs=fs, fpass=[5, 1000], tapers=[1, 1], Npairs=2000,
                      itc=1, nfft=2048)

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

        save_name = subj + '_' + condstem + '_results.mat'

        if (not os.path.isdir(respath)):
            os.mkdir(respath)
        io.savemat(respath + save_name, res)

        if not os.path.isfile(respath + save_raw_name):
            io.savemat(respath + save_raw_name, dict(x=x, subj=subj))
