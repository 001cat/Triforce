import numpy as np
from scipy import integrate

def gaussDis(mu,sig,t):
    return gaussFun(1/np.sqrt(2*np.pi*sig**2),mu,sig,t)
def gaussFun(A,mu,sig,t):
    return A*np.exp(-((t-mu)**2)/(2*sig*sig))
def gaussFit(t,f,start=[1,0,1]):
    '''Get best gaussian function fit. A,mu,sig = gaussFit(t,f)'''
    from scipy.optimize import leastsq
    def errorFun(pars,t,f):
        A,mu,sig = pars
        return gaussFun(A,mu,sig,t) - f
    A,mu,sig = leastsq(errorFun,start,args=(t,f))[0]
    return A,mu,abs(sig)

def group_into_bins(binEdges,x,y,percentile=None):
    if binEdges.ndim == 1:
        binEdgesL,binEdgesH = binEdges[:-1],binEdges[1:]
    else:
        binEdgesL,binEdgesH = binEdges[0],binEdges[1]
    binCenter = np.zeros(len(binEdgesL))
    binAvg    = np.zeros(len(binEdgesL))*np.nan
    binStd    = np.zeros(len(binEdgesL))*np.nan
    for i in range(len(binEdgesL)):
        Y = y[(x<=binEdgesH[i]) * (x>=binEdgesL[i])]
        if percentile is not None:
            Y = Y[(Y>np.percentile(Y,percentile[0])) * (Y<np.percentile(Y,percentile[1]))]
        binCenter[i] = (binEdgesH[i]+binEdgesL[i])/2
        if Y.size == 0:
            continue
        binAvg[i]    = Y.mean()
        binStd[i]    = Y.std()
    return binCenter,binAvg,binStd


def logQuad(foo,xI,xF):
    '''
    integrate in range across lots of magnitude, like 10^(-10) to 10^(10)
    '''
    xI = max(1e-30,xI)
    S = 0
    xk = np.power(10,np.arange(np.ceil(np.log10(xI)),np.log10(xF),1))
    xSeps = np.insert(np.insert(xk,-1,xF),0,xI)
    for i in range(len(xSeps)-1):
        x0,x1 = xSeps[i],xSeps[i+1]
        # S += quadpy.quad(foo,x0,x1,epsabs=1e-8,epsrel=1e-5)[0]
        S += integrate.quad(foo,x0,x1,epsabs=1e-8,epsrel=1e-5)[0]
    return S

def calCurv(x,y,t=None,N=100,debug=True):
    from scipy.interpolate import interp1d
    if t is None:
        t = x
    tNew = np.linspace(np.min(t),np.max(t),N)
    yNew = interp1d(t,y,kind='cubic')(tNew)
    xNew = interp1d(t,x,kind='cubic')(tNew)
    y_t = np.gradient(yNew,tNew)
    x_t = np.gradient(xNew,tNew)
    yy_t = np.gradient(y_t,tNew)
    xx_t = np.gradient(x_t,tNew)
    curvature_val = np.abs(xx_t * y_t - x_t * yy_t) / (x_t * x_t + y_t * y_t)**1.5
    curv = interp1d(tNew,curvature_val,kind='cubic')(t)
    if debug is True:
        print('To be finised!')
    return curv

def weighted_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    weights = weights[~np.isnan(values)]
    values  = values[~np.isnan(values)]
    if len(weights) == 0:
        return -1

    #####
    # values_Main = values[np.argmax(weights)]
    # I = abs(values - values_Main) < max(values_Main,2)
    # weights = weights[I]
    # values  = values[I]
    # if len(weights) == 0:
    #     return -1
    #####

    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return np.sqrt(variance)





