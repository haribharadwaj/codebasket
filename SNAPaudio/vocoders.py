import numpy as np
from scipy.signal import fftconvolve, hilbert, blackman, firwin, filtfilt
from scipy.signal.windows import dpss


def cams(f):
    """
    Convert a frequency (Hz) to the CAM scale.

    Parameters
    ----------
    f : float or array_like
        Frequency in Hz.

    Returns
    -------
    E : float or ndarray
        CAM scale value.
    """
    return 21.4 * np.log10(0.00437 * f + 1)


def invcams(E):
    """
    Convert a CAM scale value back to frequency in Hz.

    Parameters
    ----------
    E : float or array_like
        CAM scale value.

    Returns
    -------
    f : float or ndarray
        Frequency in Hz.
    """
    return (10 ** (E / 21.4) - 1) / 0.00437


def frames2sig(frames, tail, step):
    """
    Reconstruct a signal from its windowed frames using overlap-add.

    Parameters
    ----------
    frames : ndarray, shape (nw, numwins)
        Array containing windowed frames.
    tail : int
        Number of samples of unframed signal at the end.
    step : int
        Step size in samples between frames.

    Returns
    -------
    sig : ndarray
        Resynthesized signal.
    """
    nw, numwins = frames.shape
    sig_length = (numwins - 1) * step + tail + nw
    sig = np.zeros(sig_length)
    for k in range(numwins):
        start = k * step
        sig[start:start + nw] += frames[:, k]
    return sig


def sig2frames(sig, win, step):
    """
    Extract windowed frames from a continuous signal.

    Parameters
    ----------
    sig : array_like
        Input signal (1D array).
    win : array_like
        Window function (1D array, length nw).
    step : int
        Step size in samples.

    Returns
    -------
    frames : ndarray, shape (nw, numwins)
        Windowed frames (each column is one frame multiplied by the window).
    tail : int
        Number of samples of unframed signal at the end.
    """
    sig = np.asarray(sig)
    win = np.asarray(win)
    nw = len(win)
    numwins = (len(sig) - nw) // step + 1
    tail = len(sig) - ((numwins - 1) * step + nw)
    frames = np.zeros((nw, numwins))
    for k in range(numwins):
        start = k * step
        frames[:, k] = win * sig[start:start + nw]
    return frames, tail


def nextpow2(x):
    """
    Compute the next power of 2 greater than or equal to x.

    Parameters
    ----------
    x : float
        Input value.

    Returns
    -------
    p : int
        The smallest integer p such that 2**p >= x.
    """
    return int(np.ceil(np.log2(x)))


def gammatone_ir(cf, fs, filter_order=4, align=False, L=None):
    """
    Generate the impulse response for a gammatone filter.

    Parameters
    ----------
    cf : float
         Centre frequency in Hz.
    fs : int
         Sampling frequency in Hz.
    filter_order : int, optional
         Filter order (default: 4).
    align : bool, optional
         If True, apply phase alignment so that the impulse response peaks at t=0.
    L : int, optional
         Length of impulse response in samples. If None, it is set to 2**nextpow2(0.128*fs).

    Returns
    -------
    h : ndarray
         The gammatone impulse response.
    tc : float
         The time delay used for phase alignment (0 if align is False).
    """
    if L is None:
        L = 2 ** nextpow2(0.128 * fs)
    b = 1.019 * 24.7 * (4.37 * cf / 1000 + 1)
    tpt = (2 * np.pi) / fs
    gain = ((1.019 * b * tpt) ** filter_order) / 6
    tmp_t = np.arange(L) / fs
    tc = 0
    phase = 0
    if align:
        tc = (filter_order - 1) / (2 * np.pi * b)
        phase = -2 * np.pi * cf * tc
    h = gain * (fs ** 3) * (tmp_t ** (filter_order - 1)) * \
        np.exp(-2 * np.pi * b * tmp_t) * \
        np.cos(2 * np.pi * cf * tmp_t + phase)
    return h, tc


def gammatoneFast(x, cfs, fs=16000, align=False):
    """
    Compute basilar membrane responses using a bank of gammatone filters via FFT convolution.

    Parameters
    ----------
    x : array_like
        Input signal (1D array).
    cfs : array_like
        Centre frequencies for the filter bank (in Hz).
    fs : int, optional
        Sampling frequency (default is 16000).
    align : bool, optional
        If True, phase alignment is applied so that the impulse responses peak at t=0.

    Returns
    -------
    bm : ndarray
        Filtered signal matrix. Each column corresponds to one filter output.
    env : ndarray
        Instantaneous envelope for each filter output.
    delay : ndarray
        Delay (in samples) removed by phase alignment for each filter.
    """
    x = np.asarray(x).flatten()
    cfs = np.asarray(cfs)
    filter_order = 4
    L = 2 ** nextpow2(0.128 * fs)
    numchans = len(cfs)
    gt = np.zeros((L, numchans))
    tc = np.zeros(numchans)
    for i in range(numchans):
        h, tc_i = gammatone_ir(cfs[i], fs, filter_order=filter_order, align=align, L=L)
        gt[:, i] = h
        tc[i] = tc_i
    # Convolve x with each filter using FFT-based convolution
    Lbm = len(x) + L - 1
    bm = np.zeros((Lbm, numchans))
    for i in range(numchans):
        bm[:, i] = fftconvolve(x, gt[:, i], mode='full')
    # Compute Hilbert envelope
    env = np.abs(hilbert(bm, axis=0))
    # Delay correction: shift the outputs by the calculated delay
    delay = np.round(tc * fs).astype(int)
    for i in range(numchans):
        d = delay[i]
        bm[:, i] = np.concatenate((bm[d:, i], np.zeros(d)))
        env[:, i] = np.concatenate((env[d:, i], np.zeros(d)))
    return bm, env, delay


def rmsnormalize(x):
    """
    Normalize the input signal to have an RMS of 0.1.

    Parameters
    ----------
    x : array_like
        Input signal.

    Returns
    -------
    y : ndarray
        RMS-normalized signal.
    """
    x = np.asarray(x)
    r = np.sqrt(np.mean(x ** 2))
    return x * (0.1 / r)


def rampsound(x, fs, risetime):
    """
    Apply a ramp to a sound signal using a DPSS window.

    Parameters
    ----------
    x : array_like
        Input sound signal.
    fs : int
        Sampling frequency in Hz.
    risetime : float
        Rise time in seconds.

    Returns
    -------
    y : ndarray
        The ramped sound signal.
    """
    x = np.asarray(x)
    Nramp = int(np.ceil(fs * risetime * 2)) + 1
    w = dpss(Nramp, 1, sym=True)
    w = w - w[0]
    w = w / np.max(w)
    half = int(np.ceil(Nramp / 2))
    ones_length = max(len(x) - 2 * half, 0)
    wbig = np.concatenate((w[:half], np.ones(ones_length), w[-half:]))
    return x * wbig


def pow2db(x):
    """
    Convert power values to decibels (dB).

    Parameters
    ----------
    x : array_like
        Power values.

    Returns
    -------
    y : array_like
        Values in decibels.
    """
    return 10 * np.log10(x)


def itfs(noisy, oracle, noise, LC, fs, nfilts=128):
    """
    Apply ideal time-frequency segregation (ITFS) using binary masking on a gammatone filter bank representation.

    Parameters
    ----------
    noisy : array_like
        The target-masker mixture signal.
    oracle : array_like
        The clean target signal (same length as noisy).
    noise : array_like
        The noise signal such that oracle + noise equals the mixture.
    LC : float
        Local SNR criterion (in dB).
    fs : int
        Sampling frequency.
    nfilts : int, optional
        Number of filter bank channels (default is 128).

    Returns
    -------
    masked : ndarray
        The cleaned mixture after applying ITFS.
    """
    f_low = 80
    f_high = 8000
    cfs = invcams(np.linspace(cams(f_low), cams(f_high), nfilts))
    noisy = np.asarray(noisy).flatten()
    oracle = np.asarray(oracle).flatten()
    noise = np.asarray(noise).flatten()
    print('Extracting basilar membrane filter outputs!')
    bm_mix, _, _ = gammatoneFast(noisy, cfs, fs, align=True)
    bm_T, _, _   = gammatoneFast(oracle, cfs, fs, align=True)
    bm_N, _, _   = gammatoneFast(noise, cfs, fs, align=True)
    
    win_ms = 20
    win_overlap_ms = 10
    win = blackman(int(np.ceil(fs * win_ms * 1e-3)))
    step = int(np.ceil(fs * win_overlap_ms * 1e-3))
    
    masked = np.zeros_like(bm_mix)
    for k in range(nfilts):
        print(f'Processing filter # {k+1} / {nfilts}')
        frames_mix, _ = sig2frames(bm_mix[:, k], win, step)
        frames_T, _   = sig2frames(bm_T[:, k], win, step)
        frames_N, tail = sig2frames(bm_N[:, k], win, step)
        SNR = pow2db(np.mean(frames_T ** 2, axis=0)) - pow2db(np.mean(frames_N ** 2, axis=0))
        mask_bin = SNR > LC
        frames_masked = frames_mix * mask_bin
        masked[:, k] = frames2sig(frames_masked, tail, step)
    masked_sum = np.sum(masked, axis=1)
    return rmsnormalize(masked_sum)


def extract_envelope(signal, fs, cutoff=300, numtaps=101):
    """
    Extract the envelope of a signal by half-wave rectification followed by zero-phase low-pass filtering.

    Parameters
    ----------
    signal : array_like
        Input signal (1D array) from which to extract the envelope.
    fs : int
        Sampling frequency in Hz.
    cutoff : float, optional
        Cutoff frequency of the low-pass filter in Hz (default is 300 Hz).
    numtaps : int, optional
        Number of taps for the FIR low-pass filter (default is 101).

    Returns
    -------
    env : ndarray
        The extracted envelope.
    """
    rectified = np.maximum(signal, 0)
    taps = firwin(numtaps, cutoff, window='blackman', fs=fs)
    # Zero-phase filtering to avoid group delay
    env = filtfilt(taps, [1.0], rectified)
    return env


def noisevocode(x, fs, nfilts=24):
    """
    Apply noise vocoding using a gammatone filter bank analysis, envelope extraction,
    and resynthesis via broadband noise modulation.

    Procedure:
      1. Filter the input signal with a gammatone filter bank (channels equally spaced between 80 and 8000 Hz).
      2. For each channel, half-wave rectify the output and extract its envelope using zero-phase low-pass filtering.
      3. Generate broadband white noise and modulate it with the extracted envelope.
      4. Re-filter the modulated noise using the corresponding gammatone filter (using zero-phase filtering).
      5. RMS-match each channel to the original gammatone-filtered output.
      6. Sum across channels to produce the final noise-vocoded speech.

    Parameters
    ----------
    x : array_like
        Input sound signal (1D array).
    fs : int
        Sampling frequency in Hz.
    nfilts : int, optional
        Number of channels (default is 24). Channels are equally spaced between 80 Hz and 8000 Hz.

    Returns
    -------
    nv : ndarray
        Noise-vocoded speech signal.
    """
    f_low = 80
    f_high = 8000
    cfs = invcams(np.linspace(cams(f_low), cams(f_high), nfilts))
    
    # Analyze the input signal using the gammatone filter bank.
    bm, _, _ = gammatoneFast(x, cfs, fs, align=True)
    L = 2 ** nextpow2(0.128 * fs)  # use the same impulse response length
    
    vocoded_channels = []
    for i in range(nfilts):
        channel_signal = bm[:, i]
        # Extract envelope via half-wave rectification and zero-phase low-pass filtering.
        env = extract_envelope(channel_signal, fs, cutoff=300, numtaps=101)
        # Generate broadband white noise for this channel.
        noise_channel = np.random.randn(len(channel_signal))
        # Modulate noise with the extracted envelope.
        modulated_noise = noise_channel * env
        
        # Obtain the gammatone impulse response for this channel.
        h, _ = gammatone_ir(cfs[i], fs, filter_order=4, align=True, L=L)
        # Use zero-phase filtering (filtfilt) to re-filter the modulated noise.
        vocoded_channel = filtfilt(h, [1.0], modulated_noise)
        
        # Match the RMS of the vocoded channel to that of the original channel.
        orig_rms = np.sqrt(np.mean(channel_signal ** 2))
        vocoded_rms = np.sqrt(np.mean(vocoded_channel ** 2))
        scale = orig_rms / vocoded_rms if vocoded_rms > 0 else 0
        vocoded_channel *= scale
        
        vocoded_channels.append(vocoded_channel)
    
    # Sum across channels to produce the final noise-vocoded output.
    nv = np.sum(np.column_stack(vocoded_channels), axis=1)
    return nv
