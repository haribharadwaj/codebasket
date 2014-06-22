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

# Something wrong with I22
plvave = 0
powave = 0
for k, subj in enumerate(subjlist):
    fpath = froot + subj + '/'

    # These are so that the generated files are organized better
    respath = fpath + 'RES/'
    save_name = subj + '_plv_inducedpow_2Hz_100Hz.mat'
    print 'Loading data for subject', subj
    dat = io.loadmat(respath + save_name)

    plvave = plvave + dat['plv']
    powave = powave + db(dat['tfspec'])

plvave = plvave/len(subjlist)
powave = powave/len(subjlist)
times = dat['times'].squeeze()
freqs = dat['freqs'].squeeze()
###############################################################################
# View time-frequency plots

pl.close('all')

pl.figure()
# bline = plvave[:, :, np.logical_and(times > 0.3, times < 0.35)].mean(axis=2)
# plvave = plvave - bline[:, :, None]
pl.imshow((plvave.mean(axis=0))**2, vmin=0.01, vmax=0.035,
          extent=[times[0]-0.35, times[-1]-0.35, freqs[0], freqs[-1]],
          aspect='auto', origin='lower')
pl.xlabel('Time (s)')
pl.ylabel('Frequency (Hz)')
pl.title('Phase Locking')
pl.colorbar()
pl.show()

pl.figure()
pl.imshow(powave.mean(axis=0),
          extent=[times[0], times[-1], freqs[0], freqs[-1]],
          aspect='auto', origin='lower')
pl.xlabel('Time (s)')
pl.ylabel('Frequency (Hz)')
pl.title('Power (dB uV)')
pl.colorbar()
pl.show()

pl.figure()
plvlow = (plvave[:, freqs < 20, :]**2).mean(axis=1).T
plvlow = plvlow - (plvlow[np.logical_and(times < 0.38, times > 0.35), :].
                   mean(axis=0))
scale = np.max(plvlow[np.logical_and(times > 0, times < 0.2), :], axis=0)
plvlow = plvlow/scale

smallITD = plvlow[:, 0:2].mean(axis=1)
largeITD = plvlow[:, 3:].mean(axis=1)
pl.plot(times-0.35, smallITD, 'r', linewidth=2)
pl.plot(times-0.35, largeITD, 'b', linewidth=2)
pl.legend(('Small ITD (50, 100, 200)', 'Large ITD (400, 800)'))
pl.xlabel('Time (s)', fontsize=20)
pl.ylabel('Normalized response', fontsize=20)
pl.show()

pl.figure()
pl.plot(times - 0.35, plvlow, linewidth=2)
pl.legend(('50us', '100us', '200us', '400us', '800us'))
pl.xlabel('Time (s)', fontsize=20)
pl.ylabel('Normalized response', fontsize=20)
pl.show()
