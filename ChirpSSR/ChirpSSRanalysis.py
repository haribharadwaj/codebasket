import mne
import numpy as np
from anlffr import spectral
from scipy import io
import os
import fnmatch
from anlffr.helper import biosemi2mne as bs
from anlffr.preproc import find_blinks
from mne.preprocessing.ssp import compute_proj_epochs
from mne.time_frequency import induced_power


# Adding Files and locations
froot = '/home/hari/Documents/PythonCodes/ChirpSSR/'

subjlist = ['I41', ]

for subj in subjlist:

    fpath = froot + subj + '/'

    # These are so that the generated files are organized better
    respath = fpath + 'RES/'

    condlist = np.arange(1, 2)
    condstemlist = ['chirp20to400Hz', ]
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

                raw.info['bads'] += ['EXG3', ]
                # Filter the data for ERPs
                raw.filter(l_freq=0.5, h_freq=200, l_trans_bandwidth=0.15,
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
                raw.add_proj(blink_projs)

                # Epoching events of type
                epochs = mne.Epochs(
                    raw, eves, cond, tmin=-0.1, proj=False,
                    tmax=2.1, baseline=(-0.1, 0.0),
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

        