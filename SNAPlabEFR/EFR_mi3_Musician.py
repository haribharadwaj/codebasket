from anlffr.helper import biosemi2mne as bs
import mne
import numpy as np
from scipy import io
import os
import fnmatch
import pylab as pl

# Adding Files and locations
froot = '/home/hari/Data/Musicianship/Pilot/'

subjlist = ['p4_pu', ]

cond = [1, 2]

overwriteOld = True
for subj in subjlist:
    fpath = froot + subj + '/'

    # These are so that the generated files are organized better
    respath = fpath + 'RES/'
    condname = 'da'

    print 'Running Subject', subj, 'Condition', condname

    save_raw_name = subj + '_' + condname + '_alltrial.mat'

    if os.path.isfile(respath + save_raw_name) and not overwriteOld:
        print 'Epoched data is already available on disk!'
        print 'Loading data from:', respath + save_raw_name
        x = io.loadmat(respath + save_raw_name)['x']
        fs = 16384.0
    else:
        bdfs = fnmatch.filter(os.listdir(fpath), subj +
                              '_mi*.bdf')
        print 'No pre-epoched data found, looking for BDF files'
        print 'Viola!', len(bdfs),  'files found!'

        for k, edfname in enumerate(bdfs):
            # Load data and read event channel
            (raw, eves) = bs.importbdf(fpath + edfname, nchans=34,
                                       refchans=['EXG1', 'EXG2'],
                                       exclude=[])

            raw.info['bads'] += ['A1', 'A2', 'A30', 'A7', 'A6', 'A24',
                                 'A28', 'A29', 'A3', 'A11', 'A15', 'A16',
                                 'A17', 'A10', 'A21', 'A20', 'A25']
            raw.drop_channels(raw.info['bads'])
            # Filter the data
            raw.filter(l_freq=80., h_freq=1000., method='iir')
            # MAYBE use w/ filtering picks=np.arange(0, 17, 1))

            # raw.apply_proj()
            fs = raw.info['sfreq']

            # Epoching events of type
            epochs = mne.Epochs(
                raw, eves, cond, tmin=-0.03, proj=False,
                tmax=0.3, baseline=(-0.03, 0.0),
                reject=dict(eeg=200e-6))  # 200 regular, 50 strict

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
    y = x.mean(axis=1)
    t = epochs.times
    pl.plot(t, y[-2, :]*1e6)
    pl.xlabel('Time (s)')
    pl.ylabel('Evoked Response (uV)')
