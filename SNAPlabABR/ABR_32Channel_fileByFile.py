from anlffr.helper import biosemi2mne as bs
import mne
import numpy as np
import os
import fnmatch
import pylab as pl

# Adding Files and locations
froot = 'D:/DATA/ABR/'

subjlist = ['S107', ]

# conds = [[3, 9], [5, 10], [6, 12], [48, 144], [80, 160], [96, 192]]

conds = [[3, 9], [5, 10], [6, 12]]
for subj in subjlist:
    print 'Running Subject', subj
    for cond in conds:
        print 'Doing condition ', cond
        fpath = froot + '/' + subj + '/'
        bdfs = fnmatch.filter(os.listdir(fpath), subj + '_ABR*.bdf')

        if len(bdfs) >= 1:
            for k, bdf in enumerate(bdfs):
                edfname = fpath + bdf
                # Load data and read event channel
                extrachans = [u'GSR1', u'GSR2', u'Erg1', u'Erg2', u'Resp',
                              u'Plet', u'Temp']
                raw, eves = bs.importbdf(edfname, nchans=36,
                                         extrachans=extrachans)
                raw.set_channel_types({'EXG3': 'eeg', 'EXG4': 'eeg'})

                # Pick channels to not include in epoch rejection
                raw.info['bads'] += ['EXG3', 'EXG4', 'A1', 'A2', 'A30', 'A7',
                                     'A6', 'A24', 'A28', 'A29', 'A3', 'A11',
                                     'A15', 'A16', 'A17', 'A10', 'A21', 'A20',
                                     'A25']
                # Filter the data
                raw.filter(l_freq=70., h_freq=3000, picks=np.arange(36))

                # Epoch the data
                tmin, tmax = -0.01, 0.015
                bmin, bmax = tmin, 0.001
                epochs = mne.Epochs(raw, eves, cond, tmin=tmin, proj=False,
                                    tmax=tmax, baseline=(bmin, bmax),
                                    picks=np.arange(36),
                                    reject=dict(eeg=100e-6),
                                    verbose='WARNING')
                xtemp = epochs.get_data()
                t = epochs.times * 1e3 - 1.6  # Adjust for delay and use ms
                # Reshaping so that channels is first
                if(xtemp.shape[0] > 0):
                    xtemp = xtemp.transpose((1, 0, 2))
                    if(k == 0):
                        x = xtemp
                    else:
                        x = np.concatenate((x, xtemp), axis=1)
                else:
                    continue
        else:
            RuntimeError('No BDF files found!!')

        # Average and plot data
        goods = [28, 3, 30, 26, 4, 25, 7, 31, 22, 9, 8, 21, 11, 12, 18]
        y = x.mean(axis=1) * 1e6  # microV
        ref = y[34, :]
        z = y[goods, :].mean(axis=0) - ref
        pl.plot(t, z, linewidth=2)
        pl.xlabel('Time (ms)', fontsize=14)
        pl.ylabel('ABR (uV)', fontsize=14)
        pl.title('Left Ear', fontsize=14)
        pl.xlim((-2., 12.))
        pl.ylim((-1.0, 2.))
        ax = pl.gca()
        ax.tick_params(labelsize=14)
    # Apply legend across conditions
    pl.legend(['85 dB peSPL', '100 dB peSPL', '115 dB peSPL'])
    pl.show()
