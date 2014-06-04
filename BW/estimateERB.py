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
from scipy.integrate import quad
from scipy.interpolate import interp1d

def roexpwt(g,p,w,t):
    """roex filter as described in Oxenham and Shera (2003) equation (3)
    
    Oxenham, A. J., & Shera, C. A. (2003). Estimates of human cochlear tuning
    at low levels using forward and simultaneous masking. JARO, 4(4), 541-554.
    
    Parameters
    ----------
    g - Deviation from center frequency (relative, can be an array)
    p - slope parameter (different for upper and lower sides)
    t - Factor by which second slope is shallower than first
    w - relative weigths slopes (determines where 2nd starts to dominate)
    
    Returns
    -------
    W - filter weights (same shape as g)
    
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
    
def intRoexpwt2(g1,g2,p,w,t):
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
    
    (I,err) = quad(roexpwt,g1,g2,args = (p,w,t))
    
    return I
    
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
    - w*np.exp(-p*g2/t)*(2 + p*g2/t)/(p/t)
    
    
    lowlimit = -(1 - w)*np.exp(-p*g1)*(2 + p*g1)/p 
    - w*np.exp(-p*g1/t)*(2 + p*g1/t)/(p/t)
    
    
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
    

def fitRoexpwt(params, threshlist, lowNoiseEdge, highNoiseEdge, midear = None,
               fc = 4000, noisebw = 0.25, offFreq = True):
    """Calculates the squared error between the roex function fit and data
    
    Parameters
    ----------
    params[0]: pu - High frequency side slope parameter
    
    params[1]: pd - Low frequency side slope parameter
    
    params[2]: w - relative weigths slopes
    
    params[3]: t - Factor by which second slope is shallower than first
    
    threshlist - list of noise levels at threshold (actual data)
    
    lowNoiseEdge - Distance of near edge of noise band on low frequency side 
                   as a fraction of fc (+ve number)
                   
                   
    highNoiseEdge - Distance of near edge of noise band on high frequency side
                    as a fraction of fc (+ve number)
    
    midear - 2D array with first column being frequency and second being filter
                filter gain (dB) (if not specified, fit is done without it)
    
    fc  - Signal Frequency (default is 4 kHz), needed for middle-ear filter
    
    noisebw - Bandwidth on noise bands as a fraction of fc (default is 0.25)
    
    offFreq - Whether to allow for off-frequency listening (default True)
                   
    Returns
    -------
    sqerr - Sum of squared errors between fit and data
    
    Note
    ----
    
    Does not return the detector efficiency - To be implemented later
    
    """
    pu = params[0]
    pd = params[1]
    w = params[2]
    t = params[3]
    
    gvec = np.arange(-0.998,1.0,0.002) # 0.2% frequency steps
    
    npoints = len(threshlist) # Number of data points
    
    noise = np.zeros((npoints,gvec.shape[0]))
    
    # MAke signal vector incase needed for off-freq attenuation
    signal = np.zeros(gvec.shape)
    signal[np.argmin(abs(gvec))] = 1.0
    
    
    if offFreq:
        offsets = np.arange(-0.1,0.11,0.01)
        scale = 1.0 + offsets
        noffsets = offsets.shape[0]
        ng = gvec.shape[0]
        
        filterWts = np.zeros((noffsets,ng))
        for k, offset in enumerate(offsets):           
            paramsOffset = [pu*scale[k], pd*scale[k],w,t]
            filterWts[k,:] = getAudFilterPWT(gvec + offset, paramsOffset)
    else:
        filterWts = getAudFilterPWT(gvec, params)
    
    # If middle ear filter is specified, use it               
    if mef is not None:
        gvec_Hz = fc*(1+gvec)
        
        # Construct interpolant so that arbitrary frequencies can be queried
        getMidEarGain = interp1d(mef[:,0],mef[:,1])
        
        # Get gain in dB for our main frequency buffer
        mefGain = getMidEarGain(gvec_Hz)
        filterWts = filterWts*db2pow(mefGain)
                 
    for k, thresh in enumerate(threshlist):
        fd = lowNoiseEdge[k]
        fu = highNoiseEdge[k]      
        maskerIndex = np.logical_or(
                         np.logical_and(gvec > fu, gvec < (fu + noisebw)),
                         np.logical_and(gvec < -fd, gvec > (-fd - noisebw)))
        
        noise[k,maskerIndex] = db2pow(thresh)

    SNR = db(np.dot(signal,filterWts.T)) - db(np.dot(noise,filterWts.T))
    
    if(offFreq):
        SNRbest = SNR.max(axis = 1)
    else:
        SNRbest = SNR
    
    squerr = np.var(SNRbest)
    
    if (pu > 150 or pd > 150 or pu < 20 or pd < 20 or w > -5
        or t < 0.1 or t > 20):
            squerr +=1e3
            print 'Error = ',squerr
            return squerr
    else:
                
        print 'Error = ',squerr           
        return squerr

def getAudFilterPWT(gvec, params):
    """ Get the auditory filter weights given parameters for roex(pu,pd,w,t)
    
    Parameters
    ----------
    
    gvec - Frequency vector in units of fc 
    
    params[0]: pu - High frequency side slope parameter
    
    params[1]: pd - Low frequency side slope parameter
    
    params[2]: w - relative weigths slopes
    
    params[3]: t - Factor by which second slope is shallower than first
    
    Returns
    -------
    
    wts - Filter wts, same shape as gvec
    
    """
    
    pu = params[0]
    pd = params[1]
    w = params[2]
    t = params[3]
    
    # Making indicator functions for each side       
    lowSide = np.zeros(gvec.shape)
    lowSide[gvec < 0.0] = 1.0
    
    highSide = np.zeros(gvec.shape)
    highSide[gvec >= 0.0] = 1.0
    wts = (roexpwt(abs(gvec),pu,w,t)*highSide + 
                 roexpwt(abs(gvec),pd,w,t)*lowSide)
                 
    return wts

def db2pow(x):
    """ Converts from dB to power ratio
    
    Parameters
    ----------
    
    x - Input in dB
    
    Returns
    -------
    
    m - magnitude ratio
    """
    
    m = 10.0**(x/10.0)
    return m
    
def db(x):
    """ Converts *power* to decibels
    
    Parameters
    ----------
    
    x - Input in linear units
    
    Returns
    -------
    
    y - Equivalend in decibel units
    """
    
    y = 10*np.log10(x)
    return y


def loadData(subj, rootdir, plotOrNot = True):
    """
    Local data organization specific loading function
    """
    
    
    flist = glob(rootdir + subj + '/*_Block_??_BW_*.mat')
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
            # Ps = dat['soundLevel'].ravel()[0]
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
        
        
        if inds.sum()>0 and abs(np.diff(threshlist[inds])).max() > 4:
            thresh_unique[k] = threshlist[inds].max()
        else:
            thresh_unique[k] = threshlist[inds].mean()
        
        
    K_surrogate = thresh_unique[0]
    thresh_unique = thresh_unique - K_surrogate
    
    lowNoiseEdge = []
    highNoiseEdge = []
    nanIndex = []
    for k, asymm in enumerate(asymm_unique):
        if not np.isnan(thresh_unique[k]):
            if(asymm == 0):
                lowNoiseEdge = lowNoiseEdge + [bw_unique[k],]
                highNoiseEdge = highNoiseEdge + [bw_unique[k],]
            elif(asymm == 1.0):
                lowNoiseEdge = lowNoiseEdge + [bw_unique[k]+0.2,]
                highNoiseEdge = highNoiseEdge + [bw_unique[k],]
            elif(asymm == 2.0):
                lowNoiseEdge = lowNoiseEdge + [bw_unique[k],]
                highNoiseEdge = highNoiseEdge + [bw_unique[k]+0.2,]
        else:            
            nanIndex += [k,]
    thresh_unique = np.delete(thresh_unique,np.asarray(nanIndex))            
    return (thresh_unique, lowNoiseEdge, highNoiseEdge)

 
# Code to actually minimize

rootdir = '/home/hari/Documents/MATLAB/BW/'
mefname = '/home/hari/codebasket/BW/midear_Moore_et_al_1997.mat'

subj = 'I22'

data = loadData(subj,rootdir)


mef = io.loadmat(mefname)['midear_spline']

mef[:,1] = -mef[:,1] # Change sensitivity to gain

data = data + (mef,)

fc = 4000

minimizerOptions = dict(maxiter = 2000, disp = True,maxfev = 2000)

initialGuess = np.asarray([100,40,-30, 3.5])

# Data is taken in as a tuple of 4 items:
#  (thresholds, lowNoiseEdge, highNoiseEdge, mef)
# If the last item is left out, no middle ear filter is used for fitting


fit = minimize(fitRoexpwt,initialGuess,args = data,method = 'Nelder-Mead',
               options = minimizerOptions)
pu_best = fit['x'][0]
pd_best = fit['x'][1]
w_best = fit['x'][2]
t_best = fit['x'][3]
ERB = (intRoexpwt(0,0.4,pd_best,w_best,t_best) +
    intRoexpwt(0,0.4,pu_best,w_best,t_best))
print 'ERB = ', ERB*fc, 'Hz'   

# Plot the filter
gvec = np.arange(-1.0, 1.0, 0.002)
wts = getAudFilterPWT(gvec, fit['x'])
pl.figure()
pl.plot((1+gvec)*fc, db(wts), linewidth = 2)
pl.xlim((1000, 6000))
pl.ylim( (-60.0, 0.0))
ax = pl.gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(20)
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(20)
pl.grid(True)
pl.xlabel('Frequency (Hz)', fontsize = 20)
pl.ylabel('Filter Gain (dB)', fontsize = 20)
pl.show()