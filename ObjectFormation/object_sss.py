from mne.preprocessing.maxfilter import apply_maxfilter
import fnmatch
import os

# Adding Files and locations
# froot = '/autofs/cluster/transcend/hari/ObjectFormation/'
froot = '/autofs/cluster/transcend/MEG/objectformation/'

subjlist = ['030801', '032901', '035201', '038301', '038302', '039001',
            '042201']
paradigm = 'object'

hp_est = True


class FileConventionError(Exception):
    pass


class BadChannelListMissingError(Exception):
    pass


for subj in subjlist:
    fpath = froot + subj + '/'
    visit = fnmatch.filter(os.listdir(fpath), '?')[0]
    fpath += visit + '/'
    nruns = 3
    for run in range(1, nruns + 1):
        fifs = fnmatch.filter(os.listdir(fpath), subj + '*' + paradigm + '*' +
                              str(run) + '_raw.fif')
        print 'Viola!', len(fifs),  'files found!'
        if len(fifs) > 1:
            nruns = len(fifs)
            raise FileConventionError('ERROR: Multiple raw files '
                                      'found for same run! \n'
                                      'Check that files and file names'
                                      'stick to convention!')

        rawname = fifs[0]
        if run > 1:
            transtag = '_notrans'
        else:
            transtag = ''

        sssname = fpath + subj + '_' + paradigm + ('_' + str(run) +
                                                   transtag +
                                                   '_raw_sss.fif')
        # Maxfilter parameters
        frame = 'head'
        logname = fpath + subj + '_' + paradigm + ('_' + str(run) +
                                                   transtag +
                                                   '_maxfilter.log')
        if hp_est:
            hpname = (fpath + subj + '_' + paradigm + '_' + str(run) +
                      transtag + '_hp.txt')
            mv_hp = hpname
            mv_headpos = False  # Otherwise compensation wont be applied
        else:
            mv_hp = None
            mv_headpos = False

        mx_args = '-hpisubt amp -in 9 -out 3  -v | tee ' + logname
        badchname = fpath + 'badch.txt'
        if os.path.isfile(badchname):
            bads = open(badchname, 'r').read().strip('\n')
        else:
            raise BadChannelListMissingError('ERROR: Bad channel list '
                                             'not found!\n'
                                             'Cannot Continue!')

        # Calling maxfiler
        origin = apply_maxfilter(fpath + rawname, sssname, frame=frame,
                                 bad=bads, mv_hp=mv_hp, mv_comp='inter',
                                 mv_headpos=mv_headpos, mv_hpistep=200,
                                 mx_args=mx_args, verbose='DEBUG')
        print 'Done with file:', rawname

    # Do overflow files... doing only 1 per run for now
    for run in range(1, nruns + 1):
        fifs = fnmatch.filter(os.listdir(fpath), subj + '*' + paradigm + '*' +
                              str(run) + '_raw-1.fif')
        print 'Viola!', len(fifs),  'files found!'
        if len(fifs) > 1:
            nruns = len(fifs)
            raise FileConventionError('ERROR: Multiple overflow files '
                                      'found for same run! \n'
                                      'Check that files and file names'
                                      'stick to convention!')
        if len(fifs) == 1:
            rawname = fifs[0]
            transtag = '_notrans'
            sssname = fpath + subj + '_' + paradigm + ('_' + str(run) +
                                                       transtag +
                                                       '_raw_sss-1.fif')
            # Maxfilter parameters
            frame = 'head'
            logname = fpath + subj + '_' + paradigm + ('_' + str(run) +
                                                       transtag +
                                                       '_maxfilter-1.log')
            if hp_est:
                hpname = (fpath + subj + '_' + paradigm + '_' + str(run) +
                          transtag + '_hp-1.txt')
                mv_hp = hpname
                mv_headpos = False  # Otherwise compensation wont be applied
            else:
                mv_hp = None
                mv_headpos = False

            mv_trans = None
            mx_args = '-hpisubt amp -in 9 -out 3  -v | tee ' + logname
            badchname = fpath + 'badch.txt'
            if os.path.isfile(badchname):
                bads = open(badchname, 'r').read().strip('\n')
            else:
                raise BadChannelListMissingError('ERROR: Bad channel list '
                                                 'not found!\n'
                                                 'Cannot Continue!')

            # Calling maxfiler
            origin = apply_maxfilter(fpath + rawname, sssname, frame=frame,
                                     bad=bads, mv_hp=mv_hp, mv_comp='inter',
                                     mv_headpos=mv_headpos, mv_hpistep=200,
                                     mv_hpicons=True,
                                     mx_args=mx_args, verbose='DEBUG')
            print 'Done with file:', rawname

    # Separately run --trans to the coil definitions of run 1
    for run in range(2, nruns + 1):  # Only from run 2
        fifs = fnmatch.filter(os.listdir(fpath), subj + '*' + paradigm + '*' +
                              str(run) + '_notrans_raw_sss.fif')
        print 'Viola!', len(fifs),  'files found!'
        if len(fifs) > 1:
            nruns = len(fifs)
            raise FileConventionError('ERROR: Multiple notrans files '
                                      'found for same run! \n'
                                      'Check that files and file names'
                                      'stick to convention!')
        if len(fifs) == 1:
            sss_old = fifs[0]
            sssname = fpath + subj + '_' + paradigm + ('_' + str(run) +
                                                       '_raw_sss.fif')

            # Maxfilter parameters
            frame = 'head'
            logname = fpath + subj + '_' + paradigm + ('_' + str(run) +
                                                       '_trans_maxfilter.log')
            # Transform everything to the coil definition of run 1
            run1name = fpath + subj + '_' + paradigm + ('_1_raw.fif')
            mv_trans = run1name
            mx_args = '-force -in 9 -out 3  -v | tee ' + logname

            # Calling maxfiler
            origin = apply_maxfilter(fpath + sss_old, sssname, frame=frame,
                                     mv_trans=mv_trans, mx_args=mx_args,
                                     verbose='DEBUG')
            print 'Transformed file:', sss_old

    # Transforming overflow files... Remember doing only 1 per run for now
    for run in range(1, nruns + 1):
        fifs = fnmatch.filter(os.listdir(fpath), subj + '*' + paradigm + '*' +
                              str(run) + '_notrans_raw_sss-1.fif')
        if len(fifs) > 1:
            nruns = len(fifs)
            raise FileConventionError('ERROR: Multiple notrans overflow files'
                                      'found for same run! \n'
                                      'Check that files and file names'
                                      'stick to convention!')
        if len(fifs) == 1:
            sss_old = fifs[0]
            sssname = fpath + subj + '_' + paradigm + ('_' + str(run) +
                                                       '_raw-1.fif')
            # Maxfilter parameters
            frame = 'head'
            logname = fpath + subj + '_' + paradigm + ('_' + str(run) + '_'
                                                       'trans_maxfilter-1.log')

            # Transform everything to the coil definition of run 1
            run1name = fpath + subj + '_' + paradigm + ('_1_raw_sss.fif')
            mv_trans = run1name

            mx_args = '-force -in 9 -out 3  -v | tee ' + logname

            # Calling maxfiler
            origin = apply_maxfilter(fpath + sss_old, sssname, frame=frame,
                                     mv_trans=mv_trans, mx_args=mx_args,
                                     verbose='DEBUG')
            print 'Transformed file:', sss_old

    print '----- YAY! Done with subject ', subj
