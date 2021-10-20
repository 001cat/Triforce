import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

class GeoMap():
    # if the map is a global one, you'd better set lons like: -180,...,180
    def __init__(self,lons,lats,z):
        self.lons = lons
        self.lats = lats
        self.z    = z
        if np.any(self.lons<0):
            self.type = '-180 to 180'
        else:
            self.type = '0 to 360'
    def value(self,lon,lat):
        if (lon < 0) and self.type == '0 to 360':
            lon = lon + 360
        if (lon > 180) and self.type == '-180 to 180':
            lon = lon - 360
        if (lon-self.lons[0]) * (lon-self.lons[-1]) > 0:
            # raise ValueError('Longitude is out of range!')
            return np.nan
        if (lat-self.lats[0]) * (lat-self.lats[-1]) > 0:
            # raise ValueError('Latitude is out of range!')
            return np.nan

        i = np.where(self.lons-lon>=0)[0][0]
        j = np.where(self.lats-lat>=0)[0][0]
        z0 = self.z[j-1,i-1]
        z1 = self.z[j,i-1]
        z2 = self.z[j-1,i]
        z3 = self.z[j,i]
        Dx = self.lons[i] - self.lons[i-1]
        Dy = self.lats[j] - self.lats[j-1]
        dx = lon - self.lons[i-1]
        dy = lat - self.lats[j-1]

        z = z0+(z1-z0)*dy/Dy+(z2-z0)*dx/Dx+(z0+z3-z1-z2)*dx*dy/Dx/Dy
        return z
    def plot(self):
        XX,YY = np.meshgrid(self.lons,self.lats)
        plt.figure()
        plt.pcolormesh(XX,YY,self.z)
        

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

def savetxt(fname,*args):
    N = len(args[0])
    with open(fname,'w') as f:
        for i in range(N):
            s = ''
            for arg in args:
                s = s+f'{arg[i]} '
            s = s[:-1]+'\n'
            f.write(s)

def str2digit(s):
    if s in ['Nan','nan','None','none']:
        raise ValueError('Not supported yet')
    if s.isdigit():
        return int(s)
    try:
        float(s)
        return float(s)
    except ValueError:
        pass
    return s

def loadtxt(fname):
    data = [];init=False
    with open(fname,'r') as f:
        for line in f.readlines():
            l = line.split()
            if init is False:
                data = [ [str2digit(s)] for s in l ]
                init = True
                continue
            for s,col in zip(l,data):
                col.append(str2digit(s))
    data = [col if type(col[0]) is str else np.array(col) for col in data]
    return data



            

def get_current_memory() -> float: 
    ''' get memory usage of current process '''
    import os,psutil
    pid = os.getpid()
    p = psutil.Process(pid)
    info = p.memory_full_info()
    return info.uss / 1024. / 1024.
    # psutil.virtual_memory().percent
    # psutil.virtual_memory().available
    # psutil.swap_memory().percent

if __name__ == "__main__":
    
    # import netCDF4 as nc
    # with nc.Dataset('/home/ayu/SubdInv/models/ETOPO2v2g_f4.nc','r') as dataset:
    #     lats	= dataset.variables['y'][:].data
    #     lons	= dataset.variables['x'][:].data
    #     ages	= dataset.variables['z'][:]/100.0
    # oceAge = GeoMap(lons,lats,ages)
    # print(oceAge.value(-122,45))

    pass