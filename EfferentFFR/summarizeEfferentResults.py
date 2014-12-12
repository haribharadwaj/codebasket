import numpy as np
import pylab as pl
from scipy import io


# Adding Files and locations
froot = '/home/hari/Documents/PythonCodes/EfferentFFR/'
f331 = True
addMEM = False
if f331:
    froot = froot + 'f331/'
    f0 = 331.0
    subjlist = ['I08', 'I14', 'I33', 'I41', 'I52', 'I11', 'I13',
                'I12', 'I07']
    MEMlist = ['I02', 'I53', 'I29', 'I39', 'I54']
    if addMEM:
        subjlist += MEMlist
    ch = [3, 30, 31, 23, 22]
else:
    f0 = 100.0
    subjlist = ['I08', 'I14', 'I29', 'I33', 'I39', 'I41', 'I52', 'I11', 'I02']
    ch = [30, ]

cpca = False
if cpca:
    measureName = 'cplv'
else:
    measureName = 'S'

condstemlist = ['signalOnly', 'simultaneousNoise',
                'noise500ms_ahead', 'forwardMasking']
results = np.zeros(len(condstemlist))
normRes = np.zeros((len(subjlist), len(condstemlist) - 1))

for ks, subj in enumerate(subjlist):

    fpath = froot + subj + '/'
    # These are so that the generated files are organized better
    respath = fpath + 'RES/'

    for k, cond in enumerate(condstemlist):
        fname = respath + subj + '_' + cond + '_results.mat'
        dat = io.loadmat(fname)
        f = dat['f'].squeeze()
        if cpca:
            measure = dat[measureName].squeeze()
        else:
            measure = dat[measureName].squeeze()[ch, :].mean(axis=0)
        ind = np.argmin(np.abs(f - f0))
        results[k] = measure[ind]

    normRes[ks, 0] = np.log10(results[1]/results[0]) * 20
    normRes[ks, 1] = 0.8*np.log10(results[2]/results[1]) * 20
    normRes[ks, 2] = np.log10(results[3]/results[1]) * 20

mu = normRes.mean(axis=0)
se = normRes.std(axis=0) / len(subjlist)**0.5
pl.bar(np.arange(1, 4) - 0.25, mu, width=0.5)
pl.hold(True)
pl.errorbar(np.arange(1, 4), mu, yerr=se, fmt='ok', linewidth=3)
pl.show()
