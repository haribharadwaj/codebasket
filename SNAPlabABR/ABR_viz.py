import mne
import os
import fnmatch
import pylab as pl
from warnings import warn
import numpy as np
from anlffr.preproc import peak_finder

# Adding Files and locations
froot = 'D:/DATA/ABR/'

subjlist = ['S028', ]


def find_wave(t, y, whichwave='I'):
    """Automatically find a certain ABR wave (1, 3, or 5)

    Parameters
    ----------
    t - time array (in milliseconds)
    y - the ABR waveform (in voltage units)
    number - 0 (SP), 1 (wave I, default), 3 (wave III), or 5 (wave V)

    Returns
    -------
    locs, vals - Location of positive and negative end of peak and the values

    Note
    ----

    Currently, the waves are defined from a peak to the following trough.

    """
    tmins = dict(I=1., SP=0., III=3., V=5.)
    tmaxs = dict(I=3., SP=2., III=5., V=8.)

    tmin = tmins[whichwave]
    tmax = tmaxs[whichwave]
    i1 = np.argmin(np.abs(t - tmin))
    i2 = np.argmin(np.abs(t - tmax))

    # Peak locations and values
    plocs, pvals = peak_finder(y, (max(y) - min(y))/20.)
    valid = np.logical_and(plocs > i1, plocs < i2)
    plocs = plocs[valid].tolist()

    # Trough locations and values
    tlocs, tvals = peak_finder(y, (max(y) - min(y))/20., extrema=-1)
    valid = np.logical_and(tlocs > i1, tlocs < i2)
    tlocs = tlocs[valid].tolist()

    # Sort through peaks and troughs to get adjacent pairs (P followed by T)
    pairs = []
    while len(plocs) > 0 and len(tlocs) > 0:
        while len(tlocs) > 0 and tlocs[0] < plocs[0]:
            print len(tlocs)
            tlocs.pop(0)
        if len(tlocs) > 0:
            ploc = plocs.pop(0)
            tloc = tlocs[0]
            if len(plocs) > 0:
                ploc_next = plocs[0]
                if ploc_next > tloc:
                    tloc = tlocs.pop(0)
                    pair = [ploc, tloc]
                    pairs = pairs + [pair]
            else:
                tloc = tlocs.pop(0)
                pair = [ploc, tloc]
                pairs = pairs + [pair]

    # If there are multiple pairs return the largest peak-to-trough pair
    if len(pairs) == 0:
        warn('No peak-trough pair was found for requested wave')
        locs = []
        vals = []
        return locs, vals
    else:
        if len(pairs) == 1:
            locs = pairs[0]
            vals = [y[locs[0]], y[locs[1]]]
            return locs, vals
        else:
            warn('Multiple peak-trough pairs found. Returning the largest!')
            p2t = np.zeros(len(pairs))
            for pair in pairs:
                p2t = p2t + [y[pair[0]] - y[pair[1]], ]
            locs = pairs[np.argmax(p2t)]
            vals = [y[locs[0]], y[locs[1]]]
            return locs, vals

rebase = True
for subj in subjlist:

    fpath = froot + '/' + subj + '/'

    print 'Visualizing Subject', subj

    avefiles = fnmatch.filter(os.listdir(fpath), subj + '_ABR*-ave.fif')

    if len(avefiles) == 1:
        abrs = mne.read_evokeds(fpath + avefiles[0])
    else:
        RuntimeError('No ABR average files found!! OR too many files')

    goods = [28, 3, 30, 26, 4, 25, 7, 31, 22, 9, 8, 21, 11, 12, 18]

    # Plot data
    R = [3, 4, 5]
    L = [0, 1, 2]
    labels = ['85 dB peSPL', '100 dB peSPL', '115 dB peSPL']
    for ind, k in enumerate(L):
        abr = abrs[k]
        x = abr.data * 1e6  # microV
        t = abr.times * 1e3 - 1.6  # Adjust for delay and use milliseconds
        y = x[goods, :].mean(axis=0) - x[34, :]
        if rebase is True:
            y = y - y[(t > -2.) & (t < 0.)].mean()
        pl.plot(t, y, linewidth=2, label=labels[ind])
        locs, vals = find_wave(t, y, whichwave='I')
        # pl.hold(True)
        pl.plot(t[locs], y[locs], 'ro', markersize=6)
    pl.xlabel('Time (ms)', fontsize=14)
    pl.ylabel('ABR (uV)', fontsize=14)
    pl.title('Left Ear', fontsize=14)
    pl.xlim((-2., 13.4))
    pl.ylim((-1.0, 2.))
    ax = pl.gca()
    ax.tick_params(labelsize=14)
    pl.legend()
    pl.show()

    pl.figure()
    for ind, k in enumerate(R):
        abr = abrs[k]
        x = abr.data * 1e6  # microV
        t = abr.times * 1e3 - 1.6  # Adjust for delay and use milliseconds
        y = x[goods, :].mean(axis=0) - x[35, :]
        if rebase is True:
            y = y - y[(t > -2.) & (t < 0.)].mean()
        pl.plot(t, y, linewidth=2, label=labels[ind])
        locs, vals = find_wave(t, y, whichwave='I')
        # pl.hold(True)
        pl.plot(t[locs], y[locs], 'ro', markersize=6)
    pl.xlabel('Time (ms)', fontsize=16)
    pl.ylabel('ABR (uV)', fontsize=16)
    pl.title('Right Ear', fontsize=20)
    pl.xlim((-2., 13.4))
    pl.ylim((-1., 2.))
    ax = pl.gca()
    ax.tick_params(labelsize=14)
    pl.legend()
    pl.show()
