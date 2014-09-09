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

    condlist = np.arange(1, 16)
    f0_list = (condlist-1)*30 + 100
    condstemlist = np.array(map(str, f0_list))

    Ph_f0 = np.zeros((32, len(condlist)))

    for condind, cond in enumerate(condlist):
        condstem = condstemlist[condind]
        print 'Running Subject', subj, 'Condition', condind

        save_raw_name = subj + '_' + condstem + 'Hz_alltrial.mat'

        if os.path.isfile(respath + save_raw_name):
            print 'Epoched data is already available on disk!'
            print 'Loading data from:', respath + save_raw_name
            x = io.loadmat(respath + save_raw_name)['x']
        else:
            bdfs = fnmatch.filter(os.listdir(fpath), subj + '*FFR*.bdf')
            print 'No pre-epoched data found, looking for BDF files'
            print 'Viola!', len(bdfs),  'files found!'

            for k, edfname in enumerate(bdfs):
                # Load data and read event channel
                (raw, eves) = bs.importbdf(fpath + edfname, nchans=35,
                                           refchans=['EXG1', 'EXG2'])

                raw.info['bads'] += ['EXG3', 'A24']
                # Filter the data
                raw.filter(
                    l_freq=70, h_freq=1500, picks=np.arange(0, 34, 1))

                # raw.apply_proj()
                fs = raw.info['sfreq']

                # Epoching events of type
                epochs = mne.Epochs(
                    raw, eves, cond, tmin=-0.025, proj=False,
                    tmax=0.425, baseline=(-0.025, 0),
                    reject = dict(eeg=125e-6))

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
        fs = 4096
        params = dict(Fs=fs, fpass=[5, 1000], tapers=[1, 1], Npairs=2000,
                      itc=1)

        print 'Running Phase Estimation'
        (Ph, f) = spectral.mtphase(x, params, verbose=True)

        f0 = f0_list[condind]
        Ph_f0[:, condind] = Ph[:, np.argmin(abs(f-f0))]

        # Saving Results
        res = dict(Ph=Ph, f=f)

        save_name = subj + '_' + condstem + 'Hz_phase_results.mat'

        if (not os.path.isdir(respath)):
            os.mkdir(respath)
        io.savemat(respath + save_name, res)

        if not os.path.isfile(respath + save_raw_name):
            io.savemat(respath + save_raw_name, dict(x=x, subj=subj))

    allconds_phase_filename = subj + '_all_phase.mat'
    all_ph = dict(Ph_f0=Ph_f0, f0_list=f0_list)
    if (not os.path.isdir(respath)):
        os.mkdir(respath)
    io.savemat(respath + allconds_phase_filename, all_ph)
