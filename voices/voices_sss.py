from mne.preprocessing.maxfilter import apply_maxfilter
import fnmatch
import os

# Adding Files and locations
# froot = '/home/hari/Documents/PythonCodes/voices/'
froot = '/autofs/cluster/transcend/hari/voices/'

subjlist = ['082501', ]
paradigm = 'voices'
hp_est = True
for subj in subjlist:

    fpath = froot + subj + '/'

    fifs = fnmatch.filter(os.listdir(fpath), subj + '*raw.fif')
    print 'Viola!', len(fifs),  'files found!'
    if len(fifs) > 1:
        print 'Running files one by one...'
    for k, rawname in enumerate(fifs):
        if k > 0:
            transtag = '_notrans'
        else:
            transtag = ''

        sssname = fpath + subj + '_' + paradigm + ('_' + str(k+1) + transtag +
                                                   '_raw_sss.fif')
        # Maxfilter parameters
        frame = 'head'
        logname = fpath + subj + '_' + paradigm + ('_' + str(k+1) + transtag +
                                                   '_maxfilter.log')
        if hp_est:
            hpname = (fpath + subj + '_' + paradigm + '_' + str(k+1) +
                      transtag + '_hp.txt')
            mv_hp = hpname
            mv_headpos = True
        else:
            mv_hp = None
            mv_headpos = False

        mv_trans = None

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
                                 mv_trans=mv_trans, mx_args=mx_args,
                                 verbose='DEBUG')
        print 'Estimated head center was ', origin

    # Do overflow files... doing only 1 for now
    fifs = fnmatch.filter(os.listdir(fpath), subj + '*raw-1.fif')
    print 'Viola!', len(fifs),  'files found!'
    if len(fifs) > 1:
        nruns = len(fifs)
        print 'WARNING! Multiple raw files found'
    else:
        nruns = 1

    for k, rawname in enumerate(fifs):
        transtag = '_notrans'
        sssname = fpath + subj + '_' + paradigm + ('_' + str(k+1) + transtag +
                                                   '_raw_sss-1.fif')
        # Maxfilter parameters
        frame = 'head'
        logname = fpath + subj + '_' + paradigm + ('_' + str(k+1) + transtag +
                                                   '_maxfilter-1.log')
        if hp_est:
            hpname = (fpath + subj + '_' + paradigm + '_' + str(k+1) +
                      transtag + '_hp-1.txt')
            mv_hp = hpname
            mv_headpos = True
        else:
            mv_hp = None
            mv_headpos = False

        mv_trans = None
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
                                 mv_trans=mv_trans, mx_args=mx_args,
                                 verbose='DEBUG')
        print 'Estimated head center was ', origin


# Separately run --trans to the coil definitions of run 1
for subj in subjlist:

    fpath = froot + subj + '/'

    fifs = fnmatch.filter(os.listdir(fpath), subj + '*raw.fif')
    for k, rawname in enumerate(fifs):
        sss_old = fpath + subj + '_' + paradigm + ('_' + str(k+1) +
                                                   '_notrans_raw_sss.fif')
        sssname = fpath + subj + '_' + paradigm + ('_' + str(k+1) +
                                                   '_raw_sss.fif')

        # Transform everything to the coil definition of run 1
        if k > 0:
            run1name = fpath + subj + '_' + paradigm + ('_1_raw_sss.fif')
            mv_trans = run1name
            mx_args = '-in 9 -out 3 -v | tee ' + logname

            # Maxfilter parameters
            frame = 'head'
            logname = fpath + subj + '_' + paradigm + ('_' + str(k+1) +
                                                       '_trans_maxfilter.log')
            # Calling maxfiler
            origin = apply_maxfilter(fpath + sss_old, sssname, frame=frame,
                                     mv_trans=mv_trans, mx_args=mx_args,
                                     verbose='DEBUG')
            print 'Estimated head center was ', origin
            print '----- YAY! Done with subject ', subj

    # Do overflow files... doing only 1 for now
    fifs = fnmatch.filter(os.listdir(fpath), subj + '*raw-1.fif')
    for k, rawname in enumerate(fifs):
        sss_old = fpath + subj + '_' + paradigm + ('_' + str(k+1) +
                                                   '_notrans_raw_sss-1.fif')
        sssname = fpath + subj + '_' + paradigm + ('_' + str(k+1) +
                                                   '_raw_sss-1.fif')
        # Maxfilter parameters
        frame = 'head'
        logname = fpath + subj + '_' + paradigm + ('_' + str(k+1) +
                                                   '_trans_maxfilter-1.log')

        if k > 0:
            # Transform everything to the coil definition of run 1
            run1name = fpath + subj + '_' + paradigm + ('_1_raw_sss.fif')
            mv_trans = run1name

            mx_args = '-in 9 -out 3 -v | tee ' + logname

            # Calling maxfiler
            origin = apply_maxfilter(fpath + sss_old, sssname, frame=frame,
                                     mv_trans=mv_trans, mx_args=mx_args,
                                     verbose='DEBUG')
            print 'Estimated head center was ', origin
            print '----- YAY! Done with subject ', subj
