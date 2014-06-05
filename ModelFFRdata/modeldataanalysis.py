from anlffr.helper import biosemi2mne as bs
import mne
import numpy as np
from anlffr import spectral
from scipy import io
import os
import fnmatch

# Adding Files and locations
froot = '/home/hari/Documents/PythonCodes/ModelData/'

# List of files stems, each will be appended by run number
# Could be different for edf, fiff, eve etc.
# Use list [] and enumerate over if filenames are weird


# subjlist =['I01','I02','I03','I06','I07','I08','I09','I11','I13','I14','I15',
# 'I17_redo','I18','I19','I20','I25','I26','I27','I28','I29','I30','I37',
# 'I16','I32','I33','I34','I35','I39','I05','I36']

subjlist = ['I14']

for subj in subjlist:

    fpath = froot + subj + '/'

    # These are so that the generated files are organized better
    respath = fpath + 'RES/'

    condlist = [[1, 7], [2, 8], [3, 9], [4, 10], [5, 11], [6, 12]]
    condstemlist = ['_400Hz', '_430Hz', '_460Hz', '_490Hz', '_520Hz', '_550Hz']

    for condind, cond in enumerate(condlist):
        condstem = condstemlist[condind]
        print 'Running Subject', subj, 'Condition', condind

        save_raw_name = subj + condstem + '_alltrial.mat'

        if os.path.isfile(respath + save_raw_name):
            print 'Epoched data is already available on disk!'
            print 'Loading data from:', respath + save_raw_name
            x = io.loadmat(respath + save_raw_name)['x']
        else:
            bdfs = fnmatch.filter(os.listdir(fpath), subj + '*high*.bdf')
            print 'No pre-epoched data found, looking for BDF files'
            print 'Viola!', len(bdfs),  'files found!'

            for k, edfname in enumerate(bdfs):
                # Load data and read event channel
                (raw, eves) = bs.importbdf(fpath + edfname, nchans=35,
                                           refchans=['EXG1', 'EXG2'])

                #raw.info['bads'] += ['A14', 'A25']
                # Filter the data
                raw.filter(
                    l_freq=70, h_freq=1500, picks=np.arange(0, 34, 1))

                # raw.apply_proj()
                fs = raw.info['sfreq']

                pos = dict(up=cond[0])
                neg = dict(down=cond[1])

                # Epoching events of type 1 and 7
                epochs_pos = mne.Epochs(
                    raw, eves, pos, tmin=-0.05, proj=False,
                    tmax=0.75, baseline=(-0.05, 0),
                    reject = dict(eeg=150e-6))

                xtemp = epochs_pos.get_data()

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

                env = True  # Wheter to add negative polarities
                if env:
                    epochs_neg = mne.Epochs(
                        raw, eves, pos, tmin=-0.05, proj=False,
                        tmax=0.75, baseline=(-0.05, 0),
                        reject = dict(eeg=150e-6))
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
        fs = 4096
        params = dict(Fs=fs, fpass=[5, 1000], tapers=[1, 1], Npairs=2000,
                      itc=1)

        #        print 'Running Pairwise Spectrum Estimation'
        #       (pS,f) = spectral.mtpspec(x, params, verbose = 'DEBUG')

        print 'Running Raw Spectrum Estimation'
        (Sraw, f) = spectral.mtspecraw(x, params, verbose=True)

        print 'Running Mean Spectrum Estimation'
        (S, N, f) = spectral.mtspec(x, params, verbose=True)

        print 'Running CPCA PLV Estimation'
        (cplv, f) = spectral.mtcpca(x, params, verbose=True)

        print 'Running channel by channel PLV Estimation'
        (plv, f) = spectral.mtplv(x, params, verbose=True)

        print 'Running CPCA Power Estimation'
        (cpow, f) = spectral.mtcspec(x, params, verbose=True)

        print 'Running CPCA Power Estimation'
        (Ph, f) = spectral.mtphase(x, params, verbose=True)

        # Saving Results
        res = dict(cpow=cpow, plv=plv, cplv=cplv,
                   Sraw=Sraw, f=f, S=S, N=N, Ph=Ph)

        save_name = subj + condstem + '_results.mat'

        if (not os.path.isdir(respath)):
            os.mkdir(respath)
        io.savemat(respath + save_name, res)

        if not os.path.isfile(respath + save_raw_name):
            io.savemat(respath + save_raw_name, dict(x=x, subj=subj))
