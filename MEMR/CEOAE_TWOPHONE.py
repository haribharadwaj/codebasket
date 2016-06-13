import numpy as np
import os
import fnmatch
import pylab as pl
from scipy.io import loadmat
from anlffr.preproc import band_pass_filter, peak_finder
from anlffr.dpss import dpss_windows
from math import factorial
from statsmodels.robust.scale import stand_mad as mad
from scipy.optimize import minimize_scalar


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    """Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
       the values of the time history of the signal.
    window_size : int
       the length of the window. Must be an odd integer number.
    order : int
       the order of the polynomial used in the filtering.
       Must be less then `window_size` - 1.
    deriv: int
       the order of the derivative to compute
       (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
       the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size - 1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window,
                                                           half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs(y[1:half_window+1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve(m[::-1], y, mode='valid')


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
froot = '/Users/Hari/Dropbox/Data/MEMR/CEOAE_TWOPHONE/'

subjlist = ['I33_right_60dB']
cancelinput = True
fs = 48828.125  # Hz
input_delay = 2.2e-3  # ms
oaewin = (0., 21.)
for subj in subjlist:

    fpath = froot + subj + '/'

    print 'Running Subject', subj

    matnames = fnmatch.filter(os.listdir(fpath), subj + '*CEOAE_trial*.mat')
    print 'Hmm.. %d files found' % len(matnames)
    for kfile, matfile in enumerate(matnames):
        inputClick = loadmat(fpath + matfile)['y'].squeeze()
        Praw = loadmat(fpath + matfile)['Pcanal'].squeeze()
        P = band_pass_filter(Praw, fs, 250, 20e3, filter_length='5ms')
        # P = Praw
        locs, peaks = peak_finder(inputClick, thresh=0.2)
        win_start = np.int(oaewin[0] * fs / 1000.)
        win_end = np.int(oaewin[1] * fs / 1000.)
        win_length = win_end - win_start
        if locs.shape[0] < 117 or locs.shape[0] > 115:
            nclicks_current = locs.shape[0]
            click = np.zeros((nclicks_current, win_length))
            for k, loc in enumerate(locs):
                ind = np.arange(win_start, win_end) + loc
                click[k, :] = P[ind]
        if kfile == 0:
            clicks = click
        else:
            clicks = np.concatenate((clicks, click), axis=0)

    if cancelinput:
        # 3X Trial
        matnames = fnmatch.filter(os.listdir(fpath),
                                  subj + '*CEOAE_3xtrial*.mat')
        print 'Hmm.. %d files found' % len(matnames)
        for kfile, matfile in enumerate(matnames):
            inputClick = loadmat(fpath + matfile)['y'].squeeze()
            Praw = loadmat(fpath + matfile)['Pcanal'].squeeze()
            P = band_pass_filter(Praw, fs, 250, 20e3, filter_length='5ms')
            # P = Praw
            locs, peaks = peak_finder(inputClick, thresh=0.2, extrema=1)
            win_start = np.int(oaewin[0] * fs / 1000.)
            win_end = np.int(oaewin[1] * fs / 1000.)
            win_length = win_end - win_start
            if locs.shape[0] < 117 or locs.shape[0] > 115:
                nclicks_current = locs.shape[0]
                click3x = np.zeros((nclicks_current, win_length))
                for k, loc in enumerate(locs):
                    ind = np.arange(win_start, win_end) + loc
                    click3x[k, :] = P[ind]
            if kfile == 0:
                clicks3x = click3x
            else:
                clicks3x = np.concatenate((clicks3x, click3x), axis=0)

        # 1X + 3X trial i.e., "both" trial
        matnames = fnmatch.filter(os.listdir(fpath),
                                  subj + '*CEOAE_bothtrial*.mat')
        print 'Hmm.. %d files found' % len(matnames)
        for kfile, matfile in enumerate(matnames):
            inputClick = loadmat(fpath + matfile)['y'].squeeze()
            Praw = loadmat(fpath + matfile)['Pcanal'].squeeze()
            P = band_pass_filter(Praw, fs, 250, 20e3, filter_length='5ms')
            # P = Praw
            locs, peaks = peak_finder(inputClick, thresh=0.2, extrema=1)
            win_start = np.int(oaewin[0] * fs / 1000.)
            win_end = np.int(oaewin[1] * fs / 1000.)
            win_length = win_end - win_start
            if locs.shape[0] < 117 or locs.shape[0] > 115:
                nclicks_current = locs.shape[0]
                clickboth = np.zeros((nclicks_current, win_length))
                for k, loc in enumerate(locs):
                    ind = np.arange(win_start, win_end) + loc
                    clickboth[k, :] = P[ind]
            if kfile == 0:
                clicksboth = clickboth
            else:
                clicksboth = np.concatenate((clicksboth, clickboth), axis=0)

goods = rejecttrials(clicks)

if cancelinput:
    goods3x = rejecttrials(clicks3x)
    goodsboth = rejecttrials(clicksboth)

    goods = np.intersect1d(np.intersect1d(goods, goods3x), goodsboth)

    clicks3x_good = clicks3x[goods, :]
    ceoae3x = clicks3x_good.mean(axis=0)

    clicksboth_good = clicksboth[goods, :]
    ceoaeboth = clicksboth_good.mean(axis=0)

    clicks_good = clicks[goods, :]
    ceoae_orig = np.mean(clicks_good, axis=0).squeeze()

    t0 = oaewin[0]
    t = np.arange(0, ceoae_orig.shape[0] / fs, 1. / fs) * 1000. - t0
    oaewin2 = (4., 21.)
    # Factor should be approx 1.0
    ind = (t < 5.0) & (t > 0.)

    def match(x):
        return ((ceoae_orig + ceoae3x - ceoaeboth*x)[ind] ** 2.).sum()
    optimal = minimize_scalar(match, bounds=(0.8, 1.2))
    factor = optimal['x']
    ceoae = ceoae_orig + ceoae3x - factor*ceoaeboth

    win_start = np.int(oaewin2[0] * fs / 1000.)
    win_end = np.int(oaewin2[1] * fs / 1000.)
    ceoae = ceoae[win_start:win_end]
    t0 = t[win_start]
    t = t[win_start:win_end]
else:
    clicks_good = clicks[goods, :]
    ceoae_orig = np.mean(clicks_good, axis=0).squeeze()
    t0 = oaewin[0]
    t = np.arange(0, ceoae_orig.shape[0] / fs, 1. / fs) * 1000. - t0
    oaewin2 = (0., 21.)
    win_start = np.int(oaewin2[0] * fs / 1000.)
    win_end = np.int(oaewin2[1] * fs / 1000.)
    ceoae = ceoae_orig[win_start:win_end]
    t = t[win_start:win_end]

clicks_noise = clicks_good
clicks_noise[::2, ] *= -1.0
noise = np.mean(clicks_noise, axis=0).squeeze()
noise = noise[win_start:win_end]
Nfft = np.int(2 ** np.ceil(np.log2(ceoae.shape[0]) + 1))
f = np.arange(Nfft) * fs / Nfft
fmin, fmax = 300, 6000
w, e = dpss_windows(ceoae.shape[0], 1., 1.)
w_orig, e = dpss_windows(ceoae_orig.shape[0], 1., 1.)
S_orig = np.fft.fft(w_orig * ceoae_orig, n=Nfft).squeeze()
S = np.fft.fft(w * ceoae, n=Nfft).squeeze()
N = np.fft.fft(w * noise, n=Nfft).squeeze()
ind = np.logical_and(f >= fmin, f <= fmax)
f = f[ind]
S = S[ind]
N = N[ind]
S_orig = S_orig[ind]
f_kHz = f / 1e3

# plot absolute
pl.figure()
ax1 = pl.subplot(131)
pl.plot(f_kHz, np.log10(np.abs(S)) * 20 - 6*f_kHz, 'k', linewidth=2)
pl.hold(True)
pl.plot(f_kHz, np.log10(np.abs(N)) * 20 - 4*f_kHz, 'r--', linewidth=2)
pl.ylabel('CEOAE Magnitude (dB)', fontsize=20)
pl.xlabel('Frequency (kHz)', fontsize=20)
ax1.tick_params(labelsize=16)

ax2 = pl.subplot(132, sharex=ax1)
phase_correction = np.exp(-2 * np.pi * f * (t0 * 1e-3 - input_delay))
phi = np.unwrap(np.angle(S * phase_correction))
phi -= phi[np.argmin(np.abs(f_kHz - 0.8))]
wlen = 2 ** np.ceil(np.log2(1./np.diff(f_kHz).mean())) + 1
phi_smooth = savitzky_golay(phi, window_size=wlen, order=3)
pl.plot(f_kHz, phi, 'g', linewidth=2)
pl.hold(True)
pl.plot(f_kHz, phi_smooth, 'k', linewidth=2)
pl.ylabel('CEOAE Phase (rad)', fontsize=20)
pl.xlabel('Frequency (kHz)', fontsize=20)
ax2.tick_params(labelsize=16)

ax3 = pl.subplot(133, sharex=ax1)
group_delay = (np.diff(phi_smooth) / np.diff(f)) * 1000. / (2 * np.pi)
pl.plot(f_kHz[1:], group_delay, 'k', linewidth=2)
pl.xlabel('Frequency (kHz)', fontsize=16)
pl.ylabel('Group Delay (ms)', fontsize=20)
pl.xlim((0.5, 5.0))
ax3.tick_params(labelsize=16)
pl.show()

pl.figure()
pl.plot(t - input_delay*1e3, ceoae, linewidth=2)
pl.xlabel('Time (ms)', fontsize=16)
pl.ylabel('CEOAE', fontsize=16)
pl.show()
