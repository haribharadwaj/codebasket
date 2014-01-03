# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 22:17:31 2013

Script to analyze BW data and get Qerb
@author: Hari Bharadwaj
"""
from scipy import io
from glob import glob
import pylab as pl
import numpy as np
from scipy.optimize import minimize



def roexpwt(g,p,w,t):
    """roex filter as described in Oxenham and Shera (2003) equation (3)
    
    Oxenham, A. J., & Shera, C. A. (2003). Estimates of human cochlear tuning
    at low levels using forward and simultaneous masking. JARO, 4(4), 541-554.
    
    Parameters
    ----------
    g - Deviation from center frequency (relative)
    p - slope parameter (different for upper and lower sides)
    t - Factor by which second slope is shallower than first
    w - relative weigths slopes (determines where 2nd starts to dominate)
    
    Returns
    -------
    W - filter weighting
    
    Note
    ----
    
    To fit only a slope parameter, set w = 0 (t is immaterial then)
    To fit a roex(p,r), set t = numpy.inf 
    
    """
    w = 10.0**(w/10.0)
    W = (1 - w)*(1 + p*g)*np.exp(-1*p*g) + w*(1 + p*g/t)*np.exp(-1*p*g/t)
    return W
    
def roexpr(g,p,r):
    """roex(p,r) filter as in Patterson et al., (1982).
    
    Parameters
    ----------
    g - Deviation from center frequency (relative)
    p - slope parameter (different for upper and lower sides)
    r - Parameter detrmining transition point/skirt
    
    Returns
    -------
    W - filter weighting
    
    Note
    ----
    
    To fit only a slope parameter, set w = 0 (t is immaterial then)
    To fit a roex(p,r), set t = numpy.inf 
    
    """
    W = (1 - r)*(1 + p*g)*np.exp(-1*p*g) + r
    return W
    
def intRoexpwt(g1,g2,p,w,t):
    """ Integral of the roexpwt filter Oxenham & Shera (2003) equation (3)
    
    Parameters
    ----------
    g1, g2 - Limits of the integral in normalized terms (eg.: g1=0.1,g2=0.35)
    p - SLope parameter
    t - Factor by which second slope is shallower than first
    w - relative weigths slopes (determines where 2nd starts to dominate)
    
    Returns
    -------
    
    I - Integral of the function
    
    """
    w = 10.0**(w/10.0)
    uplimit = -(1 - w)*np.exp(-p*g2)*(2 + p*g2)/p 
    - w*np.exp(-p*g2/t)*(2 + p*g2)/(p/t)
    
    lowlimit = -(1 - w)*np.exp(-p*g1)*(2 + p*g1)/p 
    - w*np.exp(-p*g1/t)*(2 + p*g1)/(p/t)
    
    I = uplimit - lowlimit
    return I
    
     
     
    
def intRoexpr(g1,g2,p,r):
    """Calculates the integral of the roex(p,r) function from g1 to g2
    
    See Patterson et al., (1982)
    
    Parameters
    ----------
    g1, g2 - Limits of the integral in normalized terms (eg.: g1=0.1,g2=0.35)
    p - SLope parameter
    r - parameter determining transition point
    
    Returns
    --------
    I - Integral under the curve
    
    Patterson, R. D., Nimmo-Smith, I., Weber, D. L., and Milroy, R. (1982)
    “The deterioration of hearing with age: Frequency selectivity, the critical
    ratio, the audiogram, and speech threshold,” J. Acoust. Soc. Am. 72,
    1788–1803.
    
    """
    uplimit = -(1 - r)*np.exp(-p*g2)*(2 + p*g2)/p + r*g2
    lowlimit = -(1 - r)*np.exp(-p*g1)*(2 + p*g1)/p + r*g1
    
    I = uplimit - lowlimit
    return I
    

def fitRoexpr(params,bwlist,threshlist,asymmlist):
    """Calculates the squared error between the roex function fit and data
    
    Parameters
    ----------
    params[0]: pu - High frequency side slope parameter
    params[1]: pd - Low frequency side slope parameter
    params[2]: r - relative weigths of slopes
    threshlist - list of noise levels at threshold
    bwlist - list of relative notchwidths (of the nearer edge)
    asymmlist - List of asymmetries (see noiseEnergy)
    
    Returns
    -------
    sqerr - Sum of squared errors between fit and data
    
    Note
    ----
    
    Assumes that thresh is in dB re: 0 bw notch noise level
    Also assumes that noise bands have a width of 0.25*f0 on each side
    
    """
    pu = params[0]
    pd = params[1]
    r = params[2]
    
    
    constSNR = SNRroexpr(0,pu,pd,r,0,0,0)
    
    # Whether or not to assume that filter centered at fc is always used
    off_freq = True
    squerr = 0
           
    for k, bw in enumerate(bwlist):
        
        asymm = asymmlist[k]
        thresh = threshlist[k]
        
        if(off_freq):
            argslist = (pu,pd,r,bw,asymm,thresh)
            bndsf = ((-0.1, 0.1),)
            fitf = minimize(SNRroexpr,np.zeros(1),args = argslist,
                        method = 'L-BFGS-B', bounds = bndsf)
            SNR = fitf['fun']
            
        else:
            SNR = SNRroexpr(0,pu,pd,r,bw,asymm,thresh)
        
        squerr = squerr + (SNR - constSNR)**2
                
    return squerr
        
def fitRoexpwt(params,bwlist,threshlist,asymmlist):
    """Calculates the squared error between the roex function fit and data
    
    Parameters
    ----------
    params[0]: pu - High frequency side slope parameter
    params[1]: pd - Low frequency side slope parameter
    params[2]: w - relative weigths slopes
    params[3]: t - Factor by which second slope is shallower than first
    threshlist - list of noise levels at threshold
    bwlist - list of relative notchwidths (of the nearer edge)
    asymmlist - List of asymmetries (see noiseEnergy)
    
    Returns
    -------
    sqerr - Sum of squared errors between fit and data
    
    Note
    ----
    
    Assumes that thresh is in dB re: 0 bw notch noise level
    Also assumes that noise bands have a width of 0.25*f0 on each side
    
    """
    pu = params[0]
    pd = params[1]
    w = params[2]
    t = params[3]
    
    
    constSNR = SNRroexpwt(0,pu,pd,w,t,0,0,0)
    # Whether or not to assume that filter centered at fc is always used
    off_freq = True
    squerr = 0
           
    for k, bw in enumerate(bwlist):
        
        asymm = asymmlist[k]
        thresh = threshlist[k]
           
        if(off_freq):
            argslist = (pu,pd,w,t,bw,asymm,thresh)
            bndsf = ((-0.1, 0.1),)
            fitf = minimize(SNRroexpwt,np.zeros(1),args = argslist,
                        method = 'L-BFGS-B', bounds = bndsf)
            
            SNR = fitf['fun']
            
        else:
            SNR = SNRroexpwt(0,pu,pd,w,t,bw,asymm,thresh)
        
        squerr = squerr + (SNR - constSNR)**2
                
    return squerr
        

def SNRroexpwt(f,pu,pd,w,t,bw,asymm,thresh):
    """Calculates the SNR of roex(p,w,t,p) filter (even if off-frequency)
    
    Parameters
    ----------
    
    f - center frequency relative to fc (can be negative or positive)
        eg.: If fc = 4000, 3500 Hz is specified as f=(3500-4000)/4000=-0.125
        This is the filter sitting at 3500 Hz.
    pu, pd - slope parameters (different for upper and lower sides)
    
    t - Factor by which second slope is shallower than first
    w - relative weigths slopes (determines where 2nd starts to dominate)
    bw - Notch width
    asymm - Asymmetry 0, 1 or 2
    thresh - Noise level relative to bw = 0 case in dB
    
    Returns
    -------
    
    negSNR - Negative of SNR of the filter (in dB)
    
    Notes
    -----
    
    Assumes noise bands are 0.25*fc wide on each side.
    Also assumes that noise bands lie fully on upper/lower side respectively.
    This routine is useful to select filter in off-frequency listening case.
    """
    
    if( f < 0):
        sig_attenuation = roexpwt(-f,pu,w,t)
    else:
        sig_attenuation = roexpwt(f,pd,w,t)
        
    notchband = 0.25
    if(asymm == 0):
        g1u = bw
        g2u = bw + notchband
        g1l = bw
        g2l = bw + notchband
    elif(asymm == 1):
        g1u = bw
        g2u = bw + notchband
        g1l = bw + 0.2
        g2l = bw + 0.2 + notchband
    elif(asymm == 2):
        g1u = bw + 0.2
        g2u = bw + 0.2 + notchband
        g1l = bw
        g2l = bw + notchband
        
    noisepow = intRoexpwt(g1u-f,g2u-f,pu,w,t) + intRoexpwt(g1l+f,g2l+f,pd,w,t)
    
    negSNR =  thresh + db(noisepow) - db(sig_attenuation)
    return negSNR
    
def SNRroexpr(f,pu,pd,r,bw,asymm,thresh):
    """Calculates the SNR of roex(p,r) filter (even if off-frequency)
    
    Parameters
    ----------
    
    f - center frequency relative to fc (can be negative or positive)
        eg.: If fc = 4000, 3500 Hz is specified as f=(3500-4000)/4000=-0.125
        This is the filter sitting at 3500 Hz.
    pu, pd - slope parameters (different for upper and lower sides)
    
    r - Weighting function 
    
    bw - Notch width
    asymm - Asymmetry 0, 1 or 2
    thresh - Noise level relative to bw = 0 case in dB
    
    Returns
    -------
    
    negSNR - Negative of SNR of the filter (in dB)
    
    Notes
    -----
    
    Assumes noise bands are 0.25*fc wide on each side.
    Also assumes that noise bands lie fully on upper/lower side respectively.
    This routine is useful to select filter in off-frequency listening case.
    """
    
    if( f < 0):
        sig_attenuation = roexpr(-f,pu,r)
    else:
        sig_attenuation = roexpr(f,pd,r)
        
    notchband = 0.25
    if(asymm == 0):
        g1u = bw
        g2u = bw + notchband
        g1l = bw
        g2l = bw + notchband
    elif(asymm == 1):
        g1u = bw
        g2u = bw + notchband
        g1l = bw + 0.2
        g2l = bw + 0.2 + notchband
    elif(asymm == 2):
        g1u = bw + 0.2
        g2u = bw + 0.2 + notchband
        g1l = bw
        g2l = bw + notchband
        
    noisepow = intRoexpr(g1u-f,g2u-f,pu,r) + intRoexpr(g1l+f,g2l+f,pd,r)
    
    
    negSNR =  thresh + db(noisepow) - db(sig_attenuation)
    return negSNR
     
def db2mag(x):
    """ Converts from dB to magnitute ratio
    
    Parameters
    ----------
    
    x - Input in dB
    
    Returns
    -------
    
    m - magnitude ratio
    """
    
    m = 10.0**(x/20.0)
    return m
    
def db(x):
    """ Converts to decibels
    
    Parameters
    ----------
    
    x - Input in linear units
    
    Returns
    -------
    
    y - Equivalend in decibel units
    """
    
    y = 20*np.log10(x)
    return y


# Actual Code
rootdir = '/home/hari/Documents/PythonCodes/research/BW/'
subj = 'I13'

flist = glob(rootdir + subj + '/*.mat')
bwlist = []
threshlist = []
asymmlist = []

# Convert from SPL RMS to spectral level SPL
fc = 4000
noisebw = (0.25 + 0.25)*fc
spectralLevel_correction = 10*np.log10(noisebw)


for k, fname in enumerate(flist):
    dat = io.loadmat(fname)   
    if(dat.has_key('bw')):
        bw = dat['bw'].ravel()[0]
        asymm = dat['asymm'].ravel()[0]
        thresh = dat['thresh'].ravel()[0] - spectralLevel_correction
        Ps = dat['soundLevel'].ravel()[0]
        bwlist = bwlist + [bw,]
        asymmlist = asymmlist + [asymm,]
        threshlist = threshlist + [thresh,] 
        if(asymm == 0):
            pl.plot(bw,thresh,'ok',linewidth = 3, ms = 10)
            pl.hold(True)
        elif(asymm == 1):
            pl.plot(bw,thresh,'<b',linewidth = 3, ms = 10)
            pl.hold(True)
        elif(asymm == 2):
            pl.plot(bw,thresh, '>r',linewidth = 3, ms = 10)
            pl.hold(True)
        else:
            print 'Unknown notch shape!'
            
bwlist = np.asarray(bwlist)
threshlist = np.asarray(threshlist)
asymmlist = np.asarray(asymmlist)    

pl.xlabel('Relative Notch Width',fontsize = 20)
pl.ylabel('Masker Spectrum Level at Threshold (dB SPL)',fontsize = 20)
pl.xlim((-0.02, 0.22))
ax = pl.gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(20)
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(20)
pl.show()

bw_unique = np.asarray([0, 0.1, 0.2, 0.2 ,0.2])
asymm_unique = np.asarray([0, 0, 0, 1, 2])
thresh_unique = np.zeros(5)
for k, bw in enumerate(bw_unique):
    inds = np.logical_and(bwlist == bw, asymmlist == asymm_unique[k])
    thresh_unique[k] = threshlist[inds].mean()
    
K_surrogate = thresh_unique[0]
thresh_unique = thresh_unique[1:] - K_surrogate
bw_unique = bw_unique[1:]
asymm_unique = asymm_unique[1:]


# Code to actually minimize
pwt = True
if(not pwt):
    initialGuess = np.asarray([30,30,0])
    bnds = ((1,100),(1,100),(0,0.9))
    data = (bw_unique, thresh_unique, asymm_unique)
    fit = minimize(fitRoexpr,initialGuess,args = data,
                   method = 'L-BFGS-B', bounds = bnds)
    pu_best = fit['x'][0]
    pd_best = fit['x'][1]
    r_best = fit['x'][2]
    ERB = intRoexpr(0,0.4,pu_best,r_best) + intRoexpr(0,0.4,pd_best,r_best)
    print 'ERB = ', ERB*fc, 'Hz'
else:
    initialGuess = np.asarray([100,50,-25, 3.5])
    bnds = ((50,120),(35,80),(-50,-20),(3,6))
    cons = {'type': 'ineq', 'fun': lambda x:  x[0] - x[1]}
    data = (bw_unique, thresh_unique, asymm_unique)
    fit = minimize(fitRoexpwt,initialGuess,args = data,method = 'SLSQP',
                   bounds = bnds, constraints = cons)
    
    pu_best = fit['x'][0]
    pd_best = fit['x'][1]
    w_best = fit['x'][2]
    t_best = fit['x'][3]
    ERB = (intRoexpwt(0,0.4,pd_best,w_best,t_best) +
        intRoexpwt(0,0.4,pu_best,w_best,t_best))
    print 'ERB = ', ERB*fc, 'Hz'   
