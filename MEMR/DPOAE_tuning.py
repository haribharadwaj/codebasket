import numpy as np
import os
import fnmatch
import pylab as pl
from scipy.io import loadmat
from math import factorial
from anlffr.dpss import dpss_windows


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


# Adding Files and locations
# froot = '/autofs/cluster/transcend/hari/MEMR/CEOAE_TWOPHONE/'
froot = '/Users/Hari/Dropbox/Data/MEMR/DPOAE_tuning/'

subjlist = ['I13_left']
fs = 48828.125  # Hz
f2 = 4000  # Hz
for subj in subjlist:

    fpath = froot + subj + '/'

    print 'Running Subject', subj

    matnames = fnmatch.filter(os.listdir(fpath), subj + '*DPOAE*' +
                              str(np.int(f2)) + '*.mat')
    print 'Hmm.. %d files found' % len(matnames)
    DPmags = []
    f3list = []
    for kfile, matfile in enumerate(matnames):
        print 'Running file %d of %d' % (kfile+1, len(matnames))
        dat = loadmat(fpath + matfile)
        P = dat['OAE'].squeeze()
        f1 = dat['f1'][0][0]
        f3 = dat['f3'][0][0]
        L3max = dat['L3max'][0][0]
        L3min = dat['L3min'][0][0]
        f_dp = 2*f1 - f2
        OAE = np.median(P, axis=0)
        winsize = np.int(fs * 0.5)
        step = np.int(fs*0.1)
        totallength = OAE.shape[0]
        numwins = (totallength - winsize)/step + 1
        DPmag = np.zeros(numwins)
        L3 = []
        for kwin in range(numwins):
            x = OAE[(kwin * step):(kwin*step + winsize)]
            Nfft = np.int(2 ** np.ceil(np.log2(x.shape[0])))
            f = np.arange(Nfft) * fs / Nfft
            w, e = dpss_windows(x.shape[0], 1., 1.)
            w /= w.sum()
            S = np.fft.fft(w * x, n=Nfft).squeeze()
            dp_ind = dp_ind = np.argmin(np.abs(f - f_dp))
            DPmag[kwin] = 20*np.log10(np.abs(S[dp_ind]))
            L3 += [(L3max - L3min) * np.double(kwin) / numwins, ]
        DPmags += [DPmag, ]
        f3list += [f3, ]
        pl.hold(True)
        pl.plot(L3, DPmags[kfile])
pl.show()
