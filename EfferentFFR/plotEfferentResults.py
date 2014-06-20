import pylab as pl
from scipy import io

# Adding Files and locations
froot = '/home/hari/Documents/PythonCodes/EfferentFFR/'


subj = 'I13'

fpath = froot + subj + '/'

# These are so that the generated files are organized better
respath = fpath + 'RES/'

condstemlist = ['signalOnly', 'simultaneousNoise',
                'noise500ms_ahead', 'noiseOnly']

for k, cond in enumerate(condstemlist):
    fname = respath + subj + '_' + cond + '_results.mat'
    dat = io.loadmat(fname)
    f = dat['f'].squeeze()
    cpow = dat['cpow'].squeeze()
    cplv = dat['cplv'].squeeze()
    Sraw = dat['Sraw'].squeeze()

    # Plot PLV
    pl.figure(num=1)
    pl.plot(f, cpow, linewidth=2)
    pl.hold(True)
    pl.ylabel('Response Magnitude (uV^2)', fontsize=20)
    pl.xlabel('Frequency (Hz)', fontsize=20)

    # Plot power
    pl.figure(num=2)
    pl.plot(f, cplv, linewidth=2)
    pl.hold(True)
    pl.ylabel('Phase locking value', fontsize=20)
    pl.xlabel('Frequency (Hz)', fontsize=20)

    # Plot raw power
    pl.figure(num=3)
    ch = [3, 4, 25, 26, 30, 31]
    pl.plot(f, Sraw[ch, :].mean(axis=0), linewidth=2)
    pl.hold(True)
    pl.ylabel('Raw power', fontsize=20)
    pl.xlabel('Frequency (Hz)', fontsize=20)

pl.figure(num=1)
pl.legend(condstemlist)
pl.show()

pl.figure(num=2)
pl.legend(condstemlist)
pl.show()

pl.figure(num=3)
pl.legend(condstemlist)
pl.show()
