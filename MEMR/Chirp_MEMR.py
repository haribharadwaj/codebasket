from anlffr.helper import biosemi2mne as bs
import mne
import numpy as np
from scipy import io
import os
import fnmatch

# Adding Files and locations
# froot = '/home/hari/Documents/PythonCodes/MEMR/'
# froot = '/autofs/cluster/transcend/hari/MEMR/'
froot = '/Users/Hari/Dropbox/Data/MEMR/Chirp/'

# subjlist = ['I52', 'I54', 'I13', 'I14', 'I56', 'I53', 'I55', 'I41', 'I56',
#            'I51', 'I03']
subjlist = ['I13', ]
overwriteOld = True
for subj in subjlist:
    fpath = froot + subj + '/'

    # These are so that the generated files are organized better
    respath = fpath

    conds = [1, 2]
    condstems = ['chirp_cond', 'chirp_rare']
    for k, cond in enumerate(conds):
        condstem = condstems[k]
        print 'Running Subject', subj, 'Condition', cond
        save_raw_name = subj + '_' + condstem + '_alltrial.mat'

        if os.path.isfile(respath + save_raw_name) and not overwriteOld:
            print 'Epoched data is already available on disk!'
            print 'Loading data from:', respath + save_raw_name
            x = io.loadmat(respath + save_raw_name)['x']
            fs = 4096.0
        else:
            bdfs = fnmatch.filter(os.listdir(fpath), subj +
                                  '_' + '*Chirp*.bdf')

            print 'No pre-epoched data found, looking for BDF files'
            print 'Viola!', len(bdfs),  'files found!'

            for k, edfname in enumerate(bdfs):
                # Load data and read event channel
                (raw, eves) = bs.importbdf(fpath + edfname, nchans=35,
                                           refchans=['EXG1', 'EXG2'])

                raw.info['bads'] += ['EXG3', 'A1', 'A2', 'A30', 'A7', 'A6',
                                     'A24', 'A28', 'A29', 'A3', 'A11', 'A15',
                                     'A16', 'A17', 'A10', 'A21', 'A20', 'A25']
                raw.drop_channels(raw.info['bads'])

                # Filter the data
                raw.filter(
                    l_freq=15, h_freq=1000, picks=np.arange(0, 17, 1))

                # raw.apply_proj()
                fs = raw.info['sfreq']

                # Epoching events of type
                epochs = mne.Epochs(
                    raw, eves, cond, tmin=-0.5, proj=False,
                    tmax=8.2, baseline=(-0.1, 0.0), picks=np.arange(0, 17, 1),
                    reject=dict(eeg=200e-6))

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

        if not os.path.isfile(respath + save_raw_name):
            io.savemat(respath + save_raw_name, dict(x=x, subj=subj))
