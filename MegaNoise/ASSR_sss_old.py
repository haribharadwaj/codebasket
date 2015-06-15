from mne.preprocessing.maxfilter import apply_maxfilter
import fnmatch
import os

# Adding Files and locations
# froot = '/home/hari/Documents/PythonCodes/voices/'
froot = '/autofs/cluster/transcend/hari/ASSRold/'

subjlist = ['077001', ]
paradigm = 'assrold'
hp_est = True
for subj in subjlist:

    fpath = froot + subj + '/'

    fifs = fnmatch.filter(os.listdir(fpath), subj + '*raw.fif')
    print 'Viola!', len(fifs),  'files found!'
    if len(fifs) > 1:
        nruns = len(fifs)
        print 'WARNING! Multiple raw files found'
    else:
        nruns = 1

    for k, rawname in enumerate(fifs):
        sssname = fpath + subj + '_' + paradigm + ('_' + str(k+1) +
                                                   '_raw_sss.fif')
        # Maxfilter parameters
        frame = 'head'
        logname = fpath + subj + '_' + paradigm + ('_' + str(k+1) +
                                                   '_maxfilter.log')
        if hp_est:
            hpname = fpath + subj + '_' + paradigm + '_' + str(k+1) + '_hp.txt'
            mv_hp = hpname
            mv_headpos = True
        else:
            mv_hp = None
            mv_headpos = False

        mx_args = '-in 9 -out 3 -v | tee ' + logname
        badchname = fpath + 'badch.txt'
        if os.path.isfile(badchname):
            bads = open(badchname, 'r').read().strip('\n')
        else:
            print 'Sorry.. bad channel list not found! skipping subject'
            continue
        # Calling maxfiler
        origin = apply_maxfilter(fpath + rawname, sssname, frame=frame,
                                 bad=bads, mv_hp=mv_hp, mv_comp='inter',
                                 mv_headpos=mv_headpos, mv_hpicons=True,
                                 mx_args=mx_args, verbose='DEBUG')
        print 'Estimated head center was ', origin
        print '----- YAY! Done with subject ', subj
