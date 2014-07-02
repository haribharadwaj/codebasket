import numpy as np
import pylab as pl
from scipy import io


def db(x):
    """ Converts *power* to decibels

    Parameters
    ----------

    x - Input in linear units

    Returns
    -------

    y - Equivalend in decibel units
    """

    y = 10*np.log10(x)
    return y


# Adding Files and locations
froot = '/home/hari/Documents/MATLAB/SSRITD/EEG/'

# List of files stems, each will be appended by run number
# Could be different for edf, fiff, eve etc.
# Use list [] and enumerate over if filenames are weird

subjlist = [
    'I02', 'I03', 'I05', 'I06', 'I08', 'I09', 'I11', 'I13', 'I14', 'I15',
    'I17', 'I18', 'I19', 'I20', 'I25', 'I26', 'I29', 'I30', 'I33', 'I35',
    'I37', 'I36', 'I39']

top = ['I37', 'I14', 'I02', 'I19', 'I26', 'I15', 'I09', 'I03', 'I39',
       'I25', 'I13', 'I33']

bottom = ['I30', 'I36', 'I17', 'I20', 'I18', 'I29', 'I35',
          'I08', 'I05', 'I06', 'I11']

select = subjlist
plvave = 0
powave = 0
for k, subj in enumerate(select):
    fpath = froot + subj + '/'

    # These are so that the generated files are organized better
    respath = fpath + 'RES/'
    save_name = subj + '_plv_inducedpow_2Hz_100Hz.mat'
    print 'Loading data for subject', subj
    dat = io.loadmat(respath + save_name)

    plvave = plvave + dat['plv']
    powave = powave + db(dat['tfspec'])

plvave = plvave/len(select)
powave = powave/len(select)
times = dat['times'].squeeze()
freqs = dat['freqs'].squeeze()
###############################################################################
# View time-frequency plots

pl.close('all')
t0 = 0.35
pl.figure()
pl.imshow((plvave.mean(axis=0)), vmin=0.1, vmax=0.20,
          extent=[times[0] - t0, times[-1] - t0, freqs[0], freqs[-1]],
          aspect='auto', origin='lower')
pl.xlabel('Time (s)')
pl.ylabel('Frequency (Hz)')
pl.title('Phase Locking')
pl.colorbar()
pl.show()

pl.figure()
powmean = powave.mean(axis=0)
bmin, bmax = -0.2, 0.0

powbline = powmean[:, np.logical_and(times < bmax, times > bmin)].mean(axis=1)
pl.imshow((powmean.T-powbline.T).T, vmin=-2.0, vmax=1.0,
          extent=[times[0] - t0, times[-1] - t0, freqs[0], freqs[-1]],
          aspect='auto', origin='lower')
pl.xlabel('Time (s)')
pl.ylabel('Frequency (Hz)')
pl.title('Power (dB re: baseline)')
pl.colorbar()
pl.show()

pl.figure()
fmin = 5
fmax = 12
smin = 0.3
smax = 0.5
bmin = 0.35
bmax = 0.38
plvlow = (plvave[:, np.logical_and(freqs < fmax, freqs > fmin), :]**2
          ).mean(axis=1).T
plvlow = plvlow - (plvlow[np.logical_and(times < bmax, times > bmin), :].
                   mean(axis=0))

smallITD = plvlow[:, :3].mean(axis=1)
scale = np.max(smallITD[np.logical_and(times > smin, times < smax)])
smallITD = smallITD/scale
largeITD = plvlow[:, 3:].mean(axis=1)
scale = np.max(largeITD[np.logical_and(times > smin, times < smax)])
largeITD = largeITD/scale
pl.plot(times - t0, smallITD, 'r', linewidth=2)
pl.plot(times - t0, largeITD, 'b', linewidth=2)
pl.legend(('Small ITD (50, 100, 200)', 'Large ITD (400, 800)'))
pl.xlabel('Time (s)')
pl.ylabel('Normalized response')
pl.show()

pl.figure()
scale = np.max(plvlow[np.logical_and(times > smin, times < smax), :], axis=0)
pl.plot(times - t0, plvlow/scale, linewidth=2)
pl.legend(('50us', '100us', '200us', '400us', '800us'))
pl.xlabel('Time (s)')
pl.ylabel('Normalized response')
pl.show()

# Storing individual data
plv = np.zeros((len(subjlist), ) + plvave.shape)
for k, subj in enumerate(select):
    fpath = froot + subj + '/'

    # These are so that the generated files are organized better
    respath = fpath + 'RES/'
    save_name = subj + '_plv_inducedpow_2Hz_100Hz.mat'
    print 'Loading data for subject', subj
    dat = io.loadmat(respath + save_name)

    plv[k, ] = dat['plv']

fmin = 4
fmax = 10
smin = 0.3
smax = 0.5
bmin = 0.35
bmax = 0.38
plv = plv[:, :, np.logical_and(freqs > fmin, freqs < fmax), :].mean(axis=2)
small = plv[:, :3, :].mean(axis=1)
#small = small - (small[:, np.logical_and(times < bmax, times > bmin)]
#                 .mean(axis=1))[:, None]
scale = np.max(small[:, np.logical_and(times < smax, times > smin)], axis=1)
small = small/scale[:, None]
wmin, wmax = 1.17, 1.42
small1 = np.max(small[:, np.logical_and(times > wmin, times < wmax)], axis=1)
small2 = small[:, np.logical_and(times > wmin, times < wmax)].mean(axis=1)

large = plv[:, 3:, :].mean(axis=1)
#large = large - (large[:, np.logical_and(times < bmax, times > bmin)]
#                 .mean(axis=1))[:, None]
scale = np.max(large[:, np.logical_and(times < smax, times > smin)], axis=1)
large = large/scale[:, None]
large1 = np.max(large[:, np.logical_and(times > wmin, times < wmax)], axis=1)
large2 = large[:, np.logical_and(times > wmin, times < wmax)].mean(axis=1)
perf = np.array([6.729,   8.563,  18.308,  19.046,  16.732,   8.266,  19.564,
                 9.457,   6.608,   7.889,  13.516,  16.509,   7.005,  13.839,
                 8.943,   7.458,  16.509,  10.344,  10.264,  16.7,  12.319,
                 6.41,   8.943])
