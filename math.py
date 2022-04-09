import numpy as np

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