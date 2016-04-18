import numpy as np
import os
import fnmatch
import pylab as pl
from scipy.io import loadmat
from anlffr.dpss import dpss_windows
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

subjlist = ['I52_left']
fs = 48828.125  # Hz
f2c = 2000.  # [1000., 2000., 4000., 8000.]
shift = f2c * np.asarray([-0.03, -0.02, 0., 0.01, 0.02, 0.03])
f2_list = f2c + shift
rat = 1.22  # f2/f1 ratio
for subj in subjlist:

    fpath = froot + subj + '/'

    print 'Running Subject', subj
    dp_mag = []
    dp_phase = []
    prim_phase = []
    for k, f2 in enumerate(f2_list):
        fstring = subj + '_DPOAE_RAW_' + str(np.int(f2)) + '*.mat'
        matname = fnmatch.filter(os.listdir(fpath), fstring)
        print 'Hmm.. %d files found: %s' % (len(matname), matname[0])
        P = loadmat(fpath + matname[0])['OAE'].squeeze()
        good = rejecttrials(P, thresh=3.)
        x = P[good].mean(axis=0)
        Nfft = np.int(2 ** np.ceil(np.log2(x.shape[0])))
        f = np.arange(Nfft) * fs / Nfft
        w, e = dpss_windows(x.shape[0], 1., 1.)
        w /= w.sum()
        S = np.fft.fft(w * x, n=Nfft).squeeze()
        f1 = f2 / rat
        f_dp = 2 * f1 - f2
        dp_ind = np.argmin(np.abs(f - f_dp))
        f1_ind = np.argmin(np.abs(f - f1))
        f2_ind = np.argmin(np.abs(f - f2))
        dp_mag += [np.abs(S[dp_ind]), ]
        phi_dp = np.angle(S[dp_ind])
        phi_prim = 2 * np.angle(S[f1_ind]) - np.angle(S[f2_ind])
        dp_phase += [phi_dp, ]
        prim_phase += [phi_prim, ]

subtractPrim = False
if subtractPrim is True:
    # phi = np.unwrap(np.asarray(dp_phase) - np.asarray(prim_phase))
    phi = np.unwrap(np.asarray(dp_phase)) - np.unwrap(np.asarray(prim_phase))
else:
    phi = np.unwrap(np.asarray(dp_phase))

pl.plot(f2_list, phi, 'o-', linewidth=2)
pl.xlabel('Frequency (kHz)', fontsize=16)
pl.ylabel('Phase (radians)', fontsize=20)
pl.show()
