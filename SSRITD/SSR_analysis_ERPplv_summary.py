from scipy import io


# Adding Files and locations
froot = '/home/hari/Documents/MATLAB/SSRITD/EEG/'

# List of files stems, each will be appended by run number
# Could be different for edf, fiff, eve etc.
# Use list [] and enumerate over if filenames are weird

subjlist = [
    'I02', 'I03', 'I05', 'I06', 'I08', 'I09', 'I11', 'I13', 'I14', 'I15',
    'I17', 'I18', 'I19', 'I20', 'I25', 'I26', 'I29', 'I30', 'I33', 'I35',
    'I37']

# Something wrong with I36
plvave = 0
for k, subj in enumerate(subjlist):
    fpath = froot + subj + '/'

    # These are so that the generated files are organized better
    respath = fpath + 'RES/'
    save_name = subj + '_plv_inducedpow.mat'
    print 'Loading data for subject', subj
    dat = io.loadmat(respath + save_name)

    plvave = plvave + dat['plv']

plvave = plvave/len(subjlist)
times = dat['times'].squeeze()
freqs = dat['freqs'].squeeze()
###############################################################################
# View time-frequency plots
import matplotlib.pyplot as plt
plt.imshow(plvave[4]**2, vmin=0.001, vmax=0.1,
           extent=[times[0], times[-1], freqs[0], freqs[-1]],
           aspect='auto', origin='lower')
plt.xlabel('Time (s)')
plt.ylabel('Frequency (Hz)')
plt.title('Phase Locking')
plt.colorbar()
plt.show()
