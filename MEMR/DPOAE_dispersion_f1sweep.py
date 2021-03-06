import numpy as np
import os
import fnmatch
import pylab as pl
from scipy.io import loadmat
from anlffr.dpss import dpss_windows
from anlffr.preproc import band_pass_filter
from statsmodels.robust.scale import stand_mad as mad


def rejecttrials(x, thresh=5.0, bipolar=True):
    """Simple function to reject trials from numpy array data

    Parameters
    ----------
    x : ndarray, shape (n_trials, n_time)
        Data as numpy array
    thresh : float, optional, default 1.5
        Threshold in number of median absolute deviations
    bipolar : boolean, optional, default True
        If odd (even) epoch is bad, also remove next even (previous odd) trial

    Returns
    -------
    list, length n_good_trials
        original array with bad epochs removed

    """

    n_trials, n_times = x.shape
    x_max = mad(x, axis=1)
    x_med = np.median(x_max)
    x_mad = mad(x_max)
    bads = []
    for k in range(n_trials):
        if np.abs(x_max[k] - x_med) > thresh * x_mad:
            bads += [k, ]
            if bipolar is True:
                if np.mod(k, 2) == 0:
                    bads += [k + 1, ]
                else:
                    bads += [k - 1, ]
        else:
            pass
    goods = np.setdiff1d(range(n_trials), np.unique(bads))
    print '%d Good trials Found' % len(goods)
    return goods


# Adding Files and locations
# froot = '/autofs/cluster/transcend/hari/MEMR/CEOAE_TWOPHONE/'
froot = '/Users/Hari/Dropbox/Data/MEMR/DPOAE_dispersion/'

subjlist = ['I13_left']
fs = 48828.125  # Hz
names = ['a', 'b', 'c', 'd']
f2s = dict(a=1000., b=2000., c=4000., d=8000.)
dfs = dict(a=24.4, b=24.4, c=48.8, d=48.8)
rat = 1.22  # f2/f1 ratio
for name in names:
    f2 = f2s[name]
    df = dfs[name]
    rat = rat
    f1c = f2 / rat
    shift = df * np.arange(-4, 5) / 2.
    f1_list = f1c + shift

    for subj in subjlist:

        fpath = froot + subj + '/f1sweep_suppression/'

        print 'Running Subject', subj
        dp_mag = []
        dp_phase = []
        prim_phase = []
        f_dp_list = []
        for k, f1 in enumerate(f1_list):
            fstring = subj + '_DPOAE_RAW_' + str(np.int(f1)) + '*.mat'
            matname = fnmatch.filter(os.listdir(fpath), fstring)
            print 'Hmm.. %d files found: %s' % (len(matname), matname[0])
            P = loadmat(fpath + matname[0])['OAE'].squeeze()
            good = rejecttrials(P, thresh=5.)
            x = P[good].mean(axis=0)
            x = band_pass_filter(x, fs, 250, 20e3, filter_length='5ms')
            Nfft = np.int(2 ** np.ceil(np.log2(x.shape[0])))
            f = np.arange(Nfft) * fs / Nfft
            w, e = dpss_windows(x.shape[0], 1., 1.)
            w /= w.sum()
            S = np.fft.fft(w * x, n=Nfft).squeeze()
            f_dp = 2 * f1 - f2
            dp_ind = np.argmin(np.abs(f - f_dp))
            f1_ind = np.argmin(np.abs(f - f1))
            f2_ind = np.argmin(np.abs(f - f2))
            dp_mag += [np.abs(S[dp_ind]), ]
            phi_dp = np.angle(S[dp_ind])
            phi_prim = 2 * np.angle(S[f1_ind]) - np.angle(S[f2_ind])
            dp_phase += [phi_dp, ]
            prim_phase += [phi_prim, ]
            f_dp_list += [f_dp, ]

    subtractPrim = False
    if subtractPrim is True:
        phi = np.unwrap(np.asarray(dp_phase) - np.asarray(prim_phase))
        # phi = np.unwrap(np.asarray(dp_phase)) - np.unwrap(np.asarray(prim_phase))
    else:
        phi = np.unwrap(np.asarray(dp_phase))

    pl.plot(f_dp_list, phi, 'o-', linewidth=2)
    pl.xlabel('Frequency (kHz)', fontsize=16)
    pl.ylabel('Phase (radians)', fontsize=20)
    pl.show()
