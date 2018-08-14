import mne
import os
import fnmatch
import pylab as pl
import numpy as np
from anlffr.preproc import peak_finder

# Adding Files and locations
froot = 'D:/DATA/ABR/'

subjlist = ['S050', ]


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

    To fit only a slope parameter, set w = 0 (t is immaterial then)
    To fit a roex(p,r), set t = numpy.inf

    """
    tmins = dict(I=1., SP=0., III=3., V=5.)
    tmaxs = dict(I=3., SP=2., III=5., V=7.)

    tmin = tmins[whichwave]
    tmax = tmaxs[whichwave]
    i1 = np.argmin(np.abs(t - tmin))
    i2 = np.argmin(np.abs(t - tmax))

    # Peak locations and values
    plocs, pvals = peak_finder(y, (max(y) - min(y))/8.)
    valid = np.logical_and(plocs > i1, plocs < i2)
    plocs = plocs[valid]
    pvals = pvals[valid]

    # Trough locations and values
    tlocs, tvals = peak_finder(y, (max(y) - min(y))/8., extrema=-1)
    valid = np.logical_and(tlocs > i1, tlocs < i2)
    tlocs = tlocs[valid]
    tvals = tvals[valid]

    # Sort through peaks and troughs to get one pair in the right order
    if plocs.shape[0] == 0 or tlocs.shape[0] == 0:
        RuntimeWarning('No valid peak-trough pairs found')
        locs = None
        vals = None
    if plocs.shape[0] + tlocs.shape[0] > 3:
        RuntimeWarning('Too many peaks and troughs found')
        locs = np.concatenate(plocs, tlocs)
        vals = np.concatenate(pvals, tvals)
    if plocs.shape[0] + tlocs.shape[0] == 2:
        if tlocs[0] > plocs[0]:
            locs = [plocs[0], tlocs[0]]
            vals = [pvals[0], tvals[0]]
        else:
            RuntimeWarning('No valid peak-trough pairs found')
            locs = None
            vals = None
    if plocs.shape[0] + tlocs.shape[0] == 3:
        if tlocs[0] > plocs[0]:
            locs = [plocs[0], tlocs[0]]
            vals = [pvals[0], tvals[0]]
        else:
            locs = [plocs[0], tlocs[1]]
            vals = [pvals[0], tvals[1]]
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

    for k in L:
        abr = abrs[k]
        x = abr.data * 1e6  # microV
        t = abr.times * 1e3 - 1.6  # Adjust for delay and use milliseconds
        y = x[goods, :].mean(axis=0) - x[34, :]
        if rebase is True:
            y = y - y[(t > -2.) & (t < 0.)].mean()
        pl.plot(t, y, linewidth=2)
        locs, vals = find_wave(t, y)
        pl.hold(True)
        pl.plot(t[locs], y[locs], 'rx', markersize=10)
    pl.xlabel('Time (ms)', fontsize=14)
    pl.ylabel('ABR (uV)', fontsize=14)
    pl.title('Left Ear', fontsize=14)
    pl.xlim((-2., 13.4))
    pl.ylim((-1.0, 2.))
    ax = pl.gca()
    ax.tick_params(labelsize=14)
    # pl.legend(['48 dB nHL', '64 dB nHL', '80 dB nHL'])
    pl.legend(['85 dB peSPL', '100 dB peSPL', '115 dB peSPL'])
    pl.show()

    pl.figure()
    for k in R:
        abr = abrs[k]
        x = abr.data * 1e6  # microV
        t = abr.times * 1e3 - 1.6  # Adjust for delay and use milliseconds
        y = x[goods, :].mean(axis=0) - x[35, :]
        if rebase is True:
            y = y - y[(t > -2.) & (t < 0.)].mean()
        pl.plot(t, y, linewidth=2)
        locs, vals = find_wave(t, y)
        pl.hold(True)
        pl.plot(t[locs], y[locs], 'rx', markersize=10)
    pl.xlabel('Time (ms)', fontsize=16)
    pl.ylabel('ABR (uV)', fontsize=16)
    pl.title('Right Ear', fontsize=20)
    pl.xlim((-2., 13.4))
    pl.ylim((-1., 2.))
    ax = pl.gca()
    ax.tick_params(labelsize=14)
    pl.legend(['85 dB peSPL', '100 dB peSPL', '115 dB peSPL'])
    pl.show()
