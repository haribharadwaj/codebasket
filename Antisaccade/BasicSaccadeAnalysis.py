from anlffr.helper import biosemi2mne as bs
import mne
import numpy as np
from scipy import io
import os
import fnmatch
import pylab as pl


# Adding Files and locations
froot = '/Users/hari/Desktop/AntisaccadeTesting/'

subj = 'BM'

fpath = froot + subj + '/'

# print 'Running Subject', subj

bdfs = fnmatch.filter(os.listdir(fpath), subj +
                      '*saccade*.bdf')
# print 'Viola!', len(bdfs),  'files found!'

bdf = bdfs[-1]  # Use only the last one

# Load data and read event channel
raw, eves = bs.importbdf(fpath + bdf, nchans=34, exclude=[],
                           refchans=['EXG2'])

# Filter the data
raw.filter(l_freq=0., h_freq=40., method='fir')

scalings = dict(eeg=150e-6)
# raw.plot(events=eves, scalings=scalings)

# Make new event number for PL, PR, AL, AR respectively
eves2 = eves.copy()
for k, trig in enumerate(eves[:, 2]):
    if trig == 1:
        eves2[k + 1, 2] += 10
    if trig == 2:
        eves2[k + 1, 2] += 20

# Epoching events of each type
picks = [0,]  # Just EXG1 needed, already references to EXG2
# Prosaccade Left (trigger 13)
epochs = mne.Epochs(
    raw, eves2, 13, tmin=-0.3, proj=False,
    tmax=2.0, baseline=(-0.3, 0.0), picks=picks,
    reject=None)
# epochs.plot(picks=[0,], scalings=scalings)
PL = epochs.get_data().squeeze() * 1.0e3

# Prosaccade Right (trigger 14)
epochs = mne.Epochs(
    raw, eves2, 14, tmin=-0.3, proj=False,
    tmax=2.0, baseline=(-0.3, 0.0), picks=picks,
    reject=None)
# epochs.plot(picks=[0,], scalings=scalings)
PR = epochs.get_data().squeeze() * -1.0e3

# Antisaccade Left (trigger 23)
epochs = mne.Epochs(
    raw, eves2, 23, tmin=-0.3, proj=False,
    tmax=2.0, baseline=(-0.3, 0.0), picks=picks,
    reject=None)
# epochs.plot(picks=[0,], scalings=scalings)
AL = epochs.get_data().squeeze() * -1.0e3

# Antisaccade Right (trigger 24)
epochs = mne.Epochs(
    raw, eves2, 24, tmin=-0.3, proj=False,
    tmax=2.0, baseline=(-0.3, 0.0), picks=picks,
    reject=None)
# epochs.plot(picks=[0,], scalings=scalings)
AR = epochs.get_data().squeeze() * 1.0e3

t = epochs.times * 1.0e3

# Combine Pros and Antis
P = np.r_[PL, PR]
A = np.r_[AL, AR]

Pm = np.median(P, axis=0).squeeze()
Am = np.median(A, axis=0).squeeze()

# Save results as matfile
mdict = dict(t=t, P=P, A=A)
resname = fpath + subj + '_saccadeResults.mat'
io.savemat(resname, mdict)

# Plot results
pl.plot(t, Pm, 'k', linewidth=2)
pl.plot(t, Am, 'r', linewidth=2)
pl.xlabel('Time (ms)')
pl.ylabel('Horizontal EOG (mV)')
pl.legend(['Prosaccade', 'Antisaccade'], loc='lower right')
pl.show()

