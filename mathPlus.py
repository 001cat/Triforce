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


def group_into_bins(binEdges,x,y,yerr=None,percentile=None):
    if np.ma.isMaskedArray(x) or np.ma.isMaskedArray(y) or np.ma.isMaskedArray(yerr):
        mask_x = False if not np.ma.isMaskedArray(x) else x.mask
        mask_y = False if not np.ma.isMaskedArray(y) else y.mask
        mask_yerr = False if not np.ma.isMaskedArray(yerr) else yerr.mask
        mask = mask_x + mask_y + mask_yerr
        x,y = x[mask==0],y[mask==0]
        yerr = None if yerr is None else yerr[mask==0]
        x = x.compressed() if np.ma.isMaskedArray(x) else x
        y = y.compressed() if np.ma.isMaskedArray(y) else y
        yerr = yerr.compressed() if np.ma.isMaskedArray(yerr) else yerr
    
    if binEdges.ndim == 1:
        binEdgesL,binEdgesH = binEdges[:-1],binEdges[1:]
    else:
        binEdgesL,binEdgesH = binEdges[0],binEdges[1]
    binCenter = np.zeros(len(binEdgesL))
    binAvg    = np.zeros(len(binEdgesL))*np.nan
    binStd    = np.zeros(len(binEdgesL))*np.nan
    for i in range(len(binEdgesL)):
        binCenter[i] = (binEdgesH[i]+binEdgesL[i])/2

        ind = (x<=binEdgesH[i]) * (x>=binEdgesL[i]) * (~np.isnan(y))
        Y = y[ind]
        if len(Y) == 0:
            continue
        if percentile is not None:
            indPercentile = (Y>np.percentile(Y,percentile[0])) * (Y<np.percentile(Y,percentile[1]))
            Y = Y[indPercentile]
            if len(Y) == 0:
                continue
        binAvg[i]    = Y.mean()
        binStd[i]    = Y.std()

        if yerr is not None:
            Yerr = yerr[ind]
            if percentile is not None:
                Yerr = Yerr[indPercentile]
            binStd[i]    = np.sqrt(Y.std()**2 + Yerr.mean()**2)
    return binCenter,binAvg,binStd

# def group_into_bins(binEdges,x,y,yerr=None,percentile=None):
    
#     if binEdges.ndim == 1:
#         binEdgesL,binEdgesH = binEdges[:-1],binEdges[1:]
#     else:
#         binEdgesL,binEdgesH = binEdges[0],binEdges[1]
#     binCenter = np.zeros(len(binEdgesL))
#     binAvg    = np.zeros(len(binEdgesL))*np.nan
#     binStd    = np.zeros(len(binEdgesL))*np.nan
#     for i in range(len(binEdgesL)):
#         binCenter[i] = (binEdgesH[i]+binEdgesL[i])/2

#         Y = y[(x<=binEdgesH[i]) * (x>=binEdgesL[i])]
#         if percentile is not None:
#             Y = Y[(Y>np.percentile(Y,percentile[0])) * (Y<np.percentile(Y,percentile[1]))]
#         if Y.size == 0:
#             continue
#         binAvg[i]    = Y.mean()
#         binStd[i]    = Y.std()

#         if yerr is not None:
#             Yerr = yerr[(x<=binEdgesH[i]) * (x>=binEdgesL[i])]
#             if percentile is not None:
#                 Yerr = Yerr[(Y>np.percentile(Y,percentile[0])) * (Y<np.percentile(Y,percentile[1]))]
#             if Yerr.size == 0:
#                 continue
#             binStd[i]    = np.sqrt(Y.std()**2 + Yerr.mean()**2)
#     return binCenter,binAvg,binStd
# def group_into_bins(binEdges,x,y,percentile=None):
#     if binEdges.ndim == 1:
#         binEdgesL,binEdgesH = binEdges[:-1],binEdges[1:]
#     else:
#         binEdgesL,binEdgesH = binEdges[0],binEdges[1]
#     binCenter = np.zeros(len(binEdgesL))
#     binAvg    = np.zeros(len(binEdgesL))*np.nan
#     binStd    = np.zeros(len(binEdgesL))*np.nan
#     for i in range(len(binEdgesL)):
#         Y = y[(x<=binEdgesH[i]) * (x>=binEdgesL[i])]
#         if percentile is not None:
#             Y = Y[(Y>np.percentile(Y,percentile[0])) * (Y<np.percentile(Y,percentile[1]))]
#         binCenter[i] = (binEdgesH[i]+binEdgesL[i])/2
#         if Y.size == 0:
#             continue
#         binAvg[i]    = Y.mean()
#         binStd[i]    = Y.std()
#     return binCenter,binAvg,binStd

def movingAvgSmooth(x,y):
    from scipy.signal import savgol_filter
    xmin,xmax = np.min(x),np.max(x)
    bins = np.zeros((2,100))
    bins[0,:] = xmin + (xmax-xmin)*np.linspace(-0.1,0.9,100)
    bins[1,:] = xmin + (xmax-xmin)*np.linspace( 0.1,1.1,100)
    newX,binAvg,_ = group_into_bins(bins,x,y)
    newY = np.interp(newX,newX[~np.isnan(binAvg)],binAvg[~np.isnan(binAvg)])
    newY = savgol_filter(newY, 21, 3)
    return newX,newY

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





