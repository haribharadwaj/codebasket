from mne.preprocessing.maxfilter import apply_maxfilter
import fnmatch
import os

# Adding Files and locations
# froot = '/home/hari/Documents/PythonCodes/voices/'
froot = '/autofs/cluster/transcend/hari/voices/'

subjlist = ['013703', ]
paradigm = 'voices'

for subj in subjlist:

    fpath = froot + subj + '/'

    fifs = fnmatch.filter(os.listdir(fpath), subj + '*raw.fif')
    print 'Viola!', len(fifs),  'files found!'
    if len(fifs) > 1:
        print 'Wait!!.. Was expecting only one file..'
        'Going to use just one'
    rawname = fifs[0]
    sssname = fpath + subj + '_' + paradigm + 'raw_sss.fif'

    # Maxfilter parameters
    frame = 'head'
    logname = fpath + subj + '_' + paradigm + '_maxfilter.log'
    hpname = fpath + subj + '_' + paradigm + '_hp.txt'
    mx_args = '-in 9 -out 3 -v > ' + logname
    badchname = fpath + 'badch.txt'
    if os.path.isfile(badchname):
        bads = open(badchname, 'r').read().strip('\n')
    else:
        print 'Sorry.. bad channel list not found! skipping subject'
        continue

    # Calling maxfiler
    origin = apply_maxfilter(rawname, sssname, frame=frame, bad=bads,
                             mv_hp=hpname, mv_comp=True, mv_hpicons=True,
                             mx_args=mx_args, verbose='DEBUG')
    print 'Estimated head center was ', origin
    print '----- YAY! Done with subject ', subj
