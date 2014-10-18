import pylab as pl
from scipy import io

# Adding Files and locations
froot = '/home/hari/Documents/PythonCodes/EfferentFFR/'
f331 = True
noisecarr = False
shorter = True
if f331:
    froot = froot + 'f331/'
    xlim = (100, 500)
    if shorter:
        condstemlist = ['signalOnly', 'simultaneousNoise',
                        'noise500ms_ahead', 'forwardMasking']
    else:
        condstemlist = ['signalOnly', 'simultaneousNoise',
                        'noise500ms_ahead', 'noiseOnly',
                        'forwardMasking']
else:
    xlim = (70, 250)
    condstemlist = ['signalOnly', 'simultaneousNoise',
                    'noise500ms_ahead', 'noiseOnly',
                    'forwardMasking']

if noisecarr and not f331:
    froot = froot + 'NoiseCarr/'
    condstemlist = ['signalOnly', 'simultaneousNoise',
                    'noise500ms_ahead', 'forwardMasking']

subj = 'I07'

fpath = froot + subj + '/'

# These are so that the generated files are organized better
respath = fpath + 'RES/'


pl.figure()
for k, cond in enumerate(condstemlist):
    fname = respath + subj + '_' + cond + '_results.mat'
    dat = io.loadmat(fname)
    f = dat['f'].squeeze()

    cpow = dat['cpow'].squeeze()
    cplv = dat['cplv'].squeeze()
    Sraw = dat['Sraw'].squeeze()

    # Plot PLV
    ax1 = pl.subplot(3, 1, 1)
    pl.plot(f, cpow, linewidth=2)
    pl.hold(True)
    pl.ylabel('Response Magnitude (uV^2)', fontsize=16)
    pl.title(' Subject ' + subj + ' Efferent FFR results')
    pl.legend(condstemlist)

    # Plot power
    ax2 = pl.subplot(3, 1, 2, sharex=ax1)
    pl.plot(f, cplv, linewidth=2)
    pl.hold(True)
    pl.ylabel('Phase locking value', fontsize=16)

    # Plot raw power
    ax3 = pl.subplot(3, 1, 3, sharex=ax1)
    ch = [3, 4, 25, 26, 30, 31]
    pl.plot(f, Sraw[ch, :].mean(axis=0), linewidth=2)
    pl.hold(True)
    pl.ylabel('Raw power', fontsize=16)
    pl.xlabel('Frequency (Hz)', fontsize=16)

pl.xlim(xlim)
pl.show()
