import tempfile,ctypes,os,sys
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from geographiclib.geodesic import Geodesic

import subprocess
try:
    result = subprocess.run(['gmt', '--version'], stdout=subprocess.PIPE)
    if result.stdout.decode().strip().split('.')[0] >= '6':
        print('Warning: default GMT version >= 6, incompatible with GeoMap.smooth(). try /usr/bin/gmt instead')
        gmtBin = '/usr/bin/gmt'
        subprocess.run([gmtBin, '--version'], stdout=subprocess.PIPE)
except:
    raise ValueError('Package GMT not found!!!')



def randString(N):
    import random,string
    ''' Return a random string '''
    return ''.join([random.choice(string.ascii_letters + string.digits) for i in range(N)])


def savetxt(fname,*args):
    N = len(args[0])
    with open(fname,'w') as f:
        for i in range(N):
            s = ''
            for arg in args:
                s = s+f'{arg[i]} '
            s = s[:-1]+'\n'
            f.write(s)
def savetxt2(fname,cols:dict,fmt=None,header=False,header_fmt=None):
    def rowNumCheck():
        rowNums = [len(k) for k in cols.values()]
        if len(set(rowNums)) != 1:
            raise ValueError('Incompatible row lengths')
        return rowNums[0]
    def getFormats(fmt):
        try:
            iter(fmt)
            return fmt
        except:
            fmt = '' if fmt is None else fmt
            return [fmt]*len(cols)
    N = rowNumCheck()
    fmts = getFormats(fmt)
    fmts_h = getFormats(header_fmt)
    with open(fname,'w') as f:
        if header:
            s = ' '.join([f'{v:{fmt}}' for v,fmt in zip(cols.keys(),fmts_h)])+'\n'; f.write(s)
        for i in range(N):
            s = ' '.join([f'{v[i]:{fmt}}' for v,fmt in zip(cols.values(),fmts)])+'\n'; f.write(s)
def _str2digit(s):
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
                data = [ [_str2digit(s)] for s in l ]
                init = True
                continue
            for s,col in zip(l,data):
                col.append(_str2digit(s))
    data = [col if type(col[0]) is str else np.array(col) for col in data]
    return data
def loadtxt2(fname,header=False,ndarray=False):
    data = {}
    with open(fname,'r') as f:
        try:
            for h in header:
                data[h] = []
            header = False
        except:
            pass
        for line in f.readlines():
            l = line.split()
            if header is True:
                for s in l:
                    data[s] = []
                header = False
                continue
            if data == {}:
                for i in range(len(l)):
                    data[f'col_{i}'] = []
            for v,s in zip(data.values(),l):
                v.append(_str2digit(s))
    if ndarray:
        for k,v in data.items():
            data[k] = np.array(v)
    return data


def get_current_memory() -> float: 
    ''' get memory usage of current process '''
    import os,psutil
    pid = os.getpid()
    p = psutil.Process(pid)
    info = p.memory_full_info()
    # psutil.virtual_memory().percent
    # psutil.virtual_memory().available
    # psutil.swap_memory().percent
    return info.uss / 1024. / 1024.


# gaussFun and gaussFit have been moved to math.py


class GeoGrid():
    R = 6371.0
    def __init__(self,lons=[],lats=[]):
        self.lons = np.array(lons)
        self.lats = np.array(lats)
        # if np.any(self.lons<0):
        #     self._lon_range = '-180 to 180'
        # else:
        #     self._lon_range = '0 to 360'

    @property
    def dlon(self):
        return self.lons[1] - self.lons[0]
    @property
    def dlat(self):
        return self.lats[1] - self.lats[0]
    @property
    def XX(self):
        return np.meshgrid(self.lons,self.lats)[0]
    @property
    def YY(self):
        return np.meshgrid(self.lons,self.lats)[1]
    @property
    def deltaXXkm(self):
        return np.radians(self.dlon)*np.cos(np.radians(self.YY))*self.R
    @property
    def deltaYYkm(self):
        return np.radians(self.dlat)*np.ones(self.YY.shape)*self.R
    @property
    def _lon_range(self):
        if np.any(self.lons>180):
            return '0 to 360'
        return '-180 to 180'

    def interpCMD(self,latIn,lonIn,zIn,tension=0.0):
        import netCDF4 as nc4
        ''' run GMT surface command to interpolate/exterpolate scatter points to whole surface
        see Smith and Wessel, 1990. '''
        with tempfile.TemporaryDirectory() as tmpdirname:
            np.savetxt(f'{tmpdirname}/fieldIn.tmp', list(zip(lonIn,latIn,zIn)), fmt='%g')
            R = f'{self.lons.min()}/{self.lons.max()}/{self.lats.min()}/{self.lats.max()}'
            os.system(f'cd {tmpdirname} && gmt gmtset MAP_FRAME_TYPE fancy && gmt surface fieldIn.tmp -T{tension} -GfieldOut.grd -I{self.dlon}/{self.dlat} -R{R}')
            try:
                with nc4.Dataset(f'{tmpdirname}/fieldOut.grd','r') as f:
                    zOut = f.variables['z'][()]
            except:
                zOut = None
        return zOut

    def interp(self,latIn,lonIn,zIn,tension=0.0):
        from pygmt.clib import Session
        import netCDF4 as nc4
        ''' similar to interpS, but avoid possible deadlock'''
        with tempfile.TemporaryDirectory() as tmpdirname:
            tmpfileIn = f'{tmpdirname}/fieldIn.tmp'
            tmpfileOut = f'{tmpdirname}/fieldOut.tmp'
            R = f'{self.lons.min()}/{self.lons.max()}/{self.lats.min()}/{self.lats.max()}'
            lib = Session()
            libgmt = ctypes.CDLL('libgmt.so')
            function = getattr(libgmt, 'GMT_Create_Session')
            function.argtypes = [ctypes.c_char_p, ctypes.c_uint, ctypes.c_uint, ctypes.c_void_p]
            function.restype  = ctypes.c_void_p
            # session = function('pygmt-session'.encode(), 2, 2, ctypes.POINTER(ctypes.c_int)())
            session = function('pygmt-session'.encode(), 2, 2, print_func_pygmt)
            lib.session_pointer = session
            np.savetxt(tmpfileIn, list(zip(lonIn,latIn,zIn)), fmt='%g')
            arg_str = f'{tmpfileIn} -G{tmpfileOut} -I{self.dlon}/{self.dlat} -R{R} -T{tension}'
            try:
                lib.call_module(module="surface", args=arg_str)
            except:
                pass
            # except:
            #     print(f'Error in surface: {arg_str}')
            lib.__exit__(None,None,None)

            try:
                with nc4.Dataset(tmpfileOut,'r') as f:
                    zOut = f.variables['z'][()]
            except:
                zOut = None
        return zOut

    def _interp_mask_(self,latIn,lonIn,dlatCheck=None,dlonCheck=None):
        dlonCheck = self.dlon*1.01 if dlonCheck is None else dlonCheck
        dlatCheck = self.dlat*1.01 if dlatCheck is None else dlatCheck
        check = np.zeros((self.XX.shape[0],self.XX.shape[1],4))
        for lat,lon in zip(latIn,lonIn):
            ind = (abs(self.XX-lon)<=dlonCheck) * (abs(self.YY-lat)<=dlatCheck)
            check[:,:,0] += (self.XX >= lon) * (self.YY >= lat) * ind
            check[:,:,1] += (self.XX < lon)  * (self.YY >= lat) * ind
            check[:,:,2] += (self.XX < lon)  * (self.YY < lat)  * ind
            check[:,:,3] += (self.XX >= lon) * (self.YY < lat)  * ind
        check = check > 0
        mask = (check[:,:,0] * check[:,:,1] * check[:,:,2] * check[:,:,3]) == 0
        return mask



        # def checkWithIn(x,y,xData,yData,dxCheck,dyCheck):
        #     s = 0
        #     for dx in [-dxCheck,dxCheck]:
        #         for dy in [-dyCheck,dyCheck]:
        #             ind = (xData <= max(x,x+dx)) * (xData >= min(x,x+dx)) * \
        #                 (yData <= max(y,y+dy)) * (yData >= min(y,y+dy))
        #             s += ind.any()
        #     return s == 4
        # dlonCheck = self.dlon*1.01 if dlonCheck is None else dlonCheck
        # dlatCheck = self.dlat*1.01 if dlatCheck is None else dlatCheck
        # mask = np.ones(self.XX.shape,dtype=bool)
        # for i in range(mask.shape[0]):
        #     for j in range(mask.shape[1]):
        #         mask[i,j] = not checkWithIn(self.lons[j],self.lats[i],lonIn,latIn,dlonCheck,dlatCheck)
        # return mask

    def gradient(self,zIn):
        zGradx,zGrady = np.zeros(self.XX.shape),np.zeros(self.YY.shape)
        deltaLon,deltaLat = self.deltaXXkm,self.deltaYYkm
        zGradx[:,1:-1] = (zIn[:,2:] - zIn[:,:-2])/(2*deltaLon[:,1:-1])
        zGrady[1:-1,:] = (zIn[2:,:] - zIn[:-2,:])/(2*deltaLat[1:-1,:])
        zGradx[:,0] = (zIn[:,1]-zIn[:,0])/deltaLon[:,0]
        zGradx[:,-1] = (zIn[:,-1]-zIn[:,-2])/deltaLon[:,-1]
        zGrady[0,:] = (zIn[1,:]-zIn[0,:])/deltaLat[0,:]
        zGrady[-1,:] = (zIn[-1,:]-zIn[-2,:])/deltaLat[-1,:]
        return zGradx,zGrady

    def laplacian(self,zIn):
        TB = 1.0        # boundary tension, see Smith and Wessel, 1990.
        n,m = zIn.shape
        deltaLon,deltaLat = self.deltaXXkm,self.deltaYYkm
        alpha = deltaLat/deltaLon
        z = np.zeros((n+2,m+2))
        z[1:-1,1:-1] = zIn
        z[1:-1,0]    = (2*(1-TB)*z[1:-1,1]-(1-TB/2)*z[1:-1,2])/(1-3/2*TB)   # extension of left,right,upper,lower boundary
        z[1:-1,-1]   = (2*(1-TB)*z[1:-1,-2]-(1-TB/2)*z[1:-1,-3])/(1-3/2*TB)
        z[0,1:-1]    = (2*(1-TB)*z[1,1:-1]-(1-TB/2)*z[2,1:-1])/(1-3/2*TB)
        z[-1,1:-1]   = (2*(1-TB)*z[-2,1:-1]-(1-TB/2)*z[-3,1:-1])/(1-3/2*TB)
        z[0,0]       = z[0,1]+z[1,0]-z[1,1]
        z[-1,0]      = z[-1,1]+z[-2,0]-z[-2,1]
        z[0,-1]      = z[0,-2]+z[1,-1]-z[1,-2]
        z[-1,-1]     = z[-2,-1]+z[-1,-2]-z[-2,-2]
        zLapla = (z[1:-1,2:]+z[1:-1,:-2]+(alpha**2)*(z[2:,1:-1]+z[:-2,1:-1])-2*(1+alpha**2)*z[1:-1,1:-1])/deltaLon/deltaLon
        return zLapla

    def _findInd(self,lon,lat):
        j = np.where(abs(self.lons-lon)<=self.dlon/4)[0][0]
        i = np.where(abs(self.lats-lat)<=self.dlat/4)[0][0]
        return i,j
    def _findInd_linear_interp(self,lon,lat):
        if (lon < 0) and self._lon_range == '0 to 360':
            lon = lon + 360
        if (lon > 180) and self._lon_range == '-180 to 180':
            lon = lon - 360
        if (lon-self.lons[0]) * (lon-self.lons[-1]) > 0:
            # raise ValueError('Longitude is out of range!')
            return None
        if (lat-self.lats[0]) * (lat-self.lats[-1]) > 0:
            # raise ValueError('Latitude is out of range!')
            return None
        try:
            j = np.where(abs(self.lons-lon)<=self.dlon/100)[0][0]
            i = np.where(abs(self.lats-lat)<=self.dlat/100)[0][0]
            return (i,j)
        except:
            j = np.where(self.lons-lon>=0)[0][0]
            i = np.where(self.lats-lat>=0)[0][0]
            Dx = self.lons[j] - self.lons[j-1]
            Dy = self.lats[i] - self.lats[i-1]
            dx = lon - self.lons[j-1]
            dy = lat - self.lats[i-1]
            return (i,j,dx,dy,Dx,Dy)
    def copy(self):
        return self.__class__(self.lons.copy(),self.lats.copy())

@ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_void_p, ctypes.c_char_p)
def print_func_pygmt(file_pointer, message):
    message = message.decode().strip()
    print(message, file=sys.stderr, flush=True)
    return 0
class GeoMap(GeoGrid):
    # if the map is a global one, you'd better set lons like: -180,...,180
    R = 6371.0
    def __init__(self,lons=[],lats=[],z=[],mask=False):
        super().__init__(lons,lats)
        # self.lons = np.array(lons)
        # self.lats = np.array(lats)
        self.z    = np.array(z)
        # if np.any(self.lons<0):
        #     self.type = '-180 to 180'
        # else:
        #     self.type = '0 to 360'
        self.mask=mask
    @property
    def zMasked(self):
        return np.ma.masked_array(self.z,mask=self.mask)
    def interpCMD(self, latIn, lonIn, zIn, tension=0):
        self.z = super().interpCMD(latIn, lonIn, zIn, tension)
    def interp(self, latIn, lonIn, zIn, tension=0):
        self.z = super().interp(latIn, lonIn, zIn, tension)
    def gradient(self):
        return super().gradient(self.z)
    def laplacian(self):
        return super().laplacian(self.z)
    def cut(self,minlon,maxlon,minlat,maxlat):
        lat_bnd, lon_bnd = [minlat,maxlat], [minlon,maxlon]
        lat_ind = np.where((self.lats > lat_bnd[0]) & (self.lats < lat_bnd[1]))
        lon_ind = np.where((self.lons > lon_bnd[0]) & (self.lons < lon_bnd[1]))
        self.z = self.z[lat_ind[0],:]
        self.z = self.z[:,lon_ind[0]]
        self.lats = self.lats[lat_ind]
        self.lons = self.lons[lon_ind]
        pass
    
    def _lon_range_change_to(self,lon_range):
        if lon_range == self._lon_range:
            return
        if lon_range == '0 to 360':
            if np.any(np.isclose(self.lons,180)) and np.any(np.isclose(self.lons,-180)):
                self.lons = self.lons[:-1]
                self.z    = self.z[:,:-1]
            lons = self.lons + 360*(self.lons<0)
            inds = np.argsort(lons)
            lons,z = lons[inds],self.z[:,inds]
            if np.any(np.isclose(self.lons,0)) and np.any(self.lons < 0) and np.any(self.lons > 0):
                lons.append(360)
                z = np.concatenate((z,z[:,0].reshape(-1,1)),axis=1)
            self.lons,self.z = lons,z
        elif lon_range == '-180 to 180':
            if np.any(np.isclose(self.lons,0)) and np.any(np.isclose(self.lons,360)):
                self.lons = self.lons[:-1]
                self.z    = self.z[:,:-1]
            lons = self.lons - 360*(self.lons>=180)
            inds = np.argsort(lons)
            lons,z = lons[inds],self.z[:,inds]
            if np.any(np.isclose(self.lons,180)) and np.any(self.lons < 180) and np.any(self.lons > 180):
                lons.append(180)
                z = np.concatenate((z,z[:,0].reshape(-1,1)),axis=1)
            self.lons,self.z = lons,z
        else:
            raise ValueError('Wrong _lon_range: should be 0 to 360 or -180 to 180')

    def value(self,lon,lat):
        indices = self._findInd_linear_interp(lon,lat)
        if indices is None:
            return np.nan
        if len(indices) == 2:
            i,j = indices
            return self.z[i,j]
        i,j,dx,dy,Dx,Dy = indices
        if len(self.z.shape) == 2:
            z0 = self.z[i-1,j-1]
            z1 = self.z[i,j-1]
            z2 = self.z[i-1,j]
            z3 = self.z[i,j]
        elif len(self.z.shape) == 3:
            z0 = self.z[i-1,j-1,:]
            z1 = self.z[i,j-1,:]
            z2 = self.z[i-1,j,:]
            z3 = self.z[i,j,:]
        else:
            raise ValueError()
        z = z0+(z1-z0)*dy/Dy+(z2-z0)*dx/Dx+(z0+z3-z1-z2)*dx*dy/Dx/Dy
        return z


        if (lon < 0) and self._lon_range == '0 to 360':
            lon = lon + 360
        if (lon > 180) and self._lon_range == '-180 to 180':
            lon = lon - 360
        if (lon-self.lons[0]) * (lon-self.lons[-1]) > 0:
            # raise ValueError('Longitude is out of range!')
            return np.nan
        if (lat-self.lats[0]) * (lat-self.lats[-1]) > 0:
            # raise ValueError('Latitude is out of range!')
            return np.nan

        try:
            i = np.where(abs(self.lons-lon)<=self.dlon/100)[0][0]
            j = np.where(abs(self.lats-lat)<=self.dlat/100)[0][0]
            return self.z[j,i]
        except:
            i = np.where(self.lons-lon>=0)[0][0]
            j = np.where(self.lats-lat>=0)[0][0]
            if len(self.z.shape) == 2:
                z0 = self.z[j-1,i-1]
                z1 = self.z[j,i-1]
                z2 = self.z[j-1,i]
                z3 = self.z[j,i]
            elif len(self.z.shape) == 3:
                z0 = self.z[j-1,i-1,:]
                z1 = self.z[j,i-1,:]
                z2 = self.z[j-1,i,:]
                z3 = self.z[j,i,:]
            else:
                raise ValueError()
            Dx = self.lons[i] - self.lons[i-1]
            Dy = self.lats[j] - self.lats[j-1]
            dx = lon - self.lons[i-1]
            dy = lat - self.lats[j-1]
            z = z0+(z1-z0)*dy/Dy+(z2-z0)*dx/Dx+(z0+z3-z1-z2)*dx*dy/Dx/Dy
            return z
    def plot(self,area=None,cmap='rainbow',**kwargs):
        '''area = [minlon,maxlon,minlat,maxlat]'''
        if self._lon_range == '0 to 360':
            raise ValueError('Error in plotting: should convert to -180 to 180 first!')
        if area is None:
            minlon = max(self.lons[0]-(self.lons[-1]-self.lons[0])/20,-180)
            maxlon = min(self.lons[-1]+(self.lons[-1]-self.lons[0])/20,180)
            minlat = max(self.lats[0]-(self.lats[-1]-self.lats[0])/20,-90)
            maxlat = min(self.lats[-1]+(self.lats[-1]-self.lats[0])/20,90)
        else:
            minlon,maxlon,minlat,maxlat = area
        # minlon,maxlon = minlon-360*(minlon>180),maxlon-360*(maxlon>180)
        # from IPython import embed; embed()
        z = self.z if self.mask is False else np.ma.masked_array(self.z,mask=self.mask)
        import cartopy.crs as ccrs
        crsPlate = ccrs.PlateCarree()
        crs = ccrs.PlateCarree()
        fig = plt.figure()
        ax = plt.axes(projection=crs)
        ax.set_extent((minlon, maxlon, minlat, maxlat))
        ax.coastlines()
        XX,YY = np.meshgrid(self.lons,self.lats)
        im = ax.pcolormesh(XX,YY,z,cmap=cmap,**kwargs)
        gl = ax.gridlines(draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False
        # ax.set_xticks(np.arange(round(minlon),round(maxlon),(maxlon-minlon)//4), crs=ccrs.PlateCarree())
        # ax.set_yticks(np.arange(round(minlat),round(maxlat),(maxlat-minlat)//4), crs=ccrs.PlateCarree())
        fig.colorbar(im)
        return im
    def smooth(self,tension=0.0, width=50.,noExtend=True):
        lons = self.lons.round(decimals=4)
        lats = self.lats.round(decimals=4)
        tmpFname = f'tmp{randString(10)}'
        XX,YY = np.meshgrid(lons,lats)
        dlon,dlat = lons[1]-lons[0],lats[1]-lats[0]
        savetxt(f'{tmpFname}.xyz',XX.flatten(),YY.flatten(),self.z.flatten())
        with open(f'{tmpFname}.bash','w+') as f:
            REG     = f'-R{lons[0]:.4f}/{lons[-1]:.4f}/{lats[0]:.4f}/{lats[-1]:.4f}'
            f.writelines(f'{gmtBin} gmtset MAP_FRAME_TYPE fancy \n')
            f.writelines(f'{gmtBin} surface {tmpFname}.xyz -T{tension} -G{tmpFname}.grd -I{dlon:.4f}/{dlat:.4f} {REG} \n')
            f.writelines(f'{gmtBin} grdfilter {tmpFname}.grd -D4 -Fg{width} -G{tmpFname}_Smooth.grd {REG} \n')
        
        # os.system(f'bash {tmpFname}.bash')
        r = subprocess.run(['bash', f'{tmpFname}.bash'], stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        if ('Input data lie exactly on a plane.' not in r.stderr.decode()) and r.stderr.decode().split():
            print(r.stderr.decode())
        
        from netCDF4 import Dataset
        with Dataset(f'{tmpFname}_Smooth.grd') as dset:
            zSmooth = dset['z'][()]
        if noExtend:
            zSmooth[np.isnan(self.z)] = np.nan
        os.system(f'rm {tmpFname}* gmt.conf gmt.history > /dev/null 2>&1')
        return self.__class__(self.lons,self.lats,zSmooth)
    def decimate(self,geoGrid:GeoGrid):
        lons,lats = geoGrid.lons.copy(),geoGrid.lats.copy()
        newZ  = np.zeros(geoGrid.XX.shape)
        for i in range(newZ.shape[0]):
            for j in range(newZ.shape[1]):
                lat,lon = lats[i],lons[j]
                newZ[i,j]  = self.value(lon,lat)
        return self.__class__(lons,lats,newZ)
    def save(self,fname):
        np.savez_compressed(fname,lons=self.lons,lats=self.lats,z=self.z,mask=self.mask)
    def load(self,fname):
        tmp = np.load(fname,allow_pickle=True)
        lons,lats,z,mask = tmp['lons'],tmp['lats'],tmp['z'],tmp['mask']
        self.__init__(lons,lats,z,mask)
    def copy(self):
        return self.__class__(self.lons.copy(),self.lats.copy(),self.z.copy())

    def profile(self,lon1,lat1,lon2,lat2,xtype='km'):
        geoDict = Geodesic.WGS84.Inverse(lat1,lon1,lat2,lon2)
        x = np.linspace(0,geoDict['s12'],301)
        z = np.ones(x.shape)*np.nan
        for i,d in enumerate(x):
            tmp = Geodesic.WGS84.Direct(lat1,lon1,geoDict['azi1'],d)
            z[i] = self.value(tmp['lon2'],tmp['lat2'])
        if xtype.lower() == 'km':
            x = x/1000
        elif xtype == 'lat':
            x = np.linspace(lat1,lat2,301)
        elif xtype == 'lon':
            x = np.linspace(lon1,lon2,301)
        else:
            raise ValueError(f'{__class__.__name__}.{__name__}: Wrong xtype: {xtype}')
        return x,z


def loadxyz(fname):
    xyz = np.loadtxt(fname)
    lats = np.sort(np.unique(xyz[:,0]))
    lons = np.sort(np.unique(xyz[:,1]))
    geoGrid = GeoGrid(lons,lats)
    geoMaps = []
    for indZ in range(2,xyz.shape[1]):
        z = np.zeros(geoGrid.XX.shape)*np.nan
        for row in xyz:
            lon,lat = row[1],row[0]
            i,j = geoGrid._findInd(lon,lat)
            z[i,j] = row[indZ]
        geoMaps.append(GeoMap(lons,lats,z))
    return geoMaps

''' deprecated '''
class GeoSphere(GeoGrid):
    def __init__(self, lons=[], lats=[]):
        super().__init__(lons, lats)
        self.R = 6371.0
    @property
    def deltaXXkm(self):
        return np.radians(self.dlon)*np.cos(np.radians(self.YY))*self.R
    @property
    def deltaYYkm(self):
        return np.radians(self.dlat)*np.ones(self.YY.shape)*self.R
    def interpCMD(self,latIn,lonIn,zIn,tension=0.0):
        import netCDF4 as nc4
        ''' run GMT surface command to interpolate/exterpolate scatter points to whole surface
        see Smith and Wessel, 1990. '''
        with tempfile.TemporaryDirectory() as tmpdirname:
            np.savetxt(f'{tmpdirname}/fieldIn.tmp', list(zip(lonIn,latIn,zIn)), fmt='%g')
            R = f'{self.lons.min()}/{self.lons.max()}/{self.lats.min()}/{self.lats.max()}'
            os.system(f'cd {tmpdirname} && {gmtBin} gmtset MAP_FRAME_TYPE fancy && {gmtBin} surface fieldIn.tmp -T{tension} -GfieldOut.grd -I{self.dlon}/{self.dlat} -R{R}')
            try:
                with nc4.Dataset(f'{tmpdirname}/fieldOut.grd','r') as f:
                    zOut = f.variables['z'][()]
            except:
                zOut = None
        return GeoMap(self.lons,self.lats,zOut)

    def interp(self,latIn,lonIn,zIn,tension=0.0):
        from pygmt.clib import Session
        import netCDF4 as nc4
        ''' similar to interpS, but avoid possible deadlock'''
        with tempfile.TemporaryDirectory() as tmpdirname:
            tmpfileIn = f'{tmpdirname}/fieldIn.tmp'
            tmpfileOut = f'{tmpdirname}/fieldOut.tmp'
            R = f'{self.lons.min()}/{self.lons.max()}/{self.lats.min()}/{self.lats.max()}'
            lib = Session()
            libgmt = ctypes.CDLL('libgmt.so')
            function = getattr(libgmt, 'GMT_Create_Session')
            function.argtypes = [ctypes.c_char_p, ctypes.c_uint, ctypes.c_uint, ctypes.c_void_p]
            function.restype  = ctypes.c_void_p
            # session = function('pygmt-session'.encode(), 2, 2, ctypes.POINTER(ctypes.c_int)())
            session = function('pygmt-session'.encode(), 2, 2, print_func_pygmt)
            lib.session_pointer = session
            np.savetxt(tmpfileIn, list(zip(lonIn,latIn,zIn)), fmt='%g')
            arg_str = f'{tmpfileIn} -G{tmpfileOut} -I{self.dlon}/{self.dlat} -R{R} -T{tension}'
            try:
                lib.call_module(module="surface", args=arg_str)
            except:
                pass
            # except:
            #     print(f'Error in surface: {arg_str}')
            lib.__exit__(None,None,None)

            try:
                with nc4.Dataset(tmpfileOut,'r') as f:
                    zOut = f.variables['z'][()]
            except:
                zOut = None
        return GeoMap(self.lons,self.lats,zOut)

    def gradient(self,zIn):
        zGradx,zGrady = np.zeros(self.XX.shape),np.zeros(self.YY.shape)
        deltaLon,deltaLat = self.deltaXXkm,self.deltaYYkm
        zGradx[:,1:-1] = (zIn[:,2:] - zIn[:,:-2])/(2*deltaLon[:,1:-1])
        zGrady[1:-1,:] = (zIn[2:,:] - zIn[:-2,:])/(2*deltaLat[1:-1,:])
        zGradx[:,0] = (zIn[:,1]-zIn[:,0])/deltaLon[:,0]
        zGradx[:,-1] = (zIn[:,-1]-zIn[:,-2])/deltaLon[:,-1]
        zGrady[0,:] = (zIn[1,:]-zIn[0,:])/deltaLat[0,:]
        zGrady[-1,:] = (zIn[-1,:]-zIn[-2,:])/deltaLat[-1,:]
        return zGradx,zGrady

    def laplacian(self,zIn):
        TB = 1.0        # boundary tension, see Smith and Wessel, 1990.
        n,m = zIn.shape
        deltaLon,deltaLat = self.deltaXXkm,self.deltaYYkm
        alpha = deltaLat/deltaLon
        z = np.zeros((n+2,m+2))
        z[1:-1,1:-1] = zIn
        z[1:-1,0]    = (2*(1-TB)*z[1:-1,1]-(1-TB/2)*z[1:-1,2])/(1-3/2*TB)   # extension of left,right,upper,lower boundary
        z[1:-1,-1]   = (2*(1-TB)*z[1:-1,-2]-(1-TB/2)*z[1:-1,-3])/(1-3/2*TB)
        z[0,1:-1]    = (2*(1-TB)*z[1,1:-1]-(1-TB/2)*z[2,1:-1])/(1-3/2*TB)
        z[-1,1:-1]   = (2*(1-TB)*z[-2,1:-1]-(1-TB/2)*z[-3,1:-1])/(1-3/2*TB)
        z[0,0]       = z[0,1]+z[1,0]-z[1,1]
        z[-1,0]      = z[-1,1]+z[-2,0]-z[-2,1]
        z[0,-1]      = z[0,-2]+z[1,-1]-z[1,-2]
        z[-1,-1]     = z[-2,-1]+z[-1,-2]-z[-2,-2]
        zLapla = (z[1:-1,2:]+z[1:-1,:-2]+(alpha**2)*(z[2:,1:-1]+z[:-2,1:-1])-2*(1+alpha**2)*z[1:-1,1:-1])/deltaLon/deltaLon
        return zLapla


if __name__ == "__main__":
    
    import netCDF4 as nc
    with nc.Dataset('/home/ayu/Packages/etopo2/ETOPO2v2g_f4.nc','r') as dataset:
        lats	= dataset.variables['y'][:].data
        lons	= dataset.variables['x'][:].data
        topo	= dataset.variables['z'][:]/1000.0
    topoMap = GeoMap(lons,lats,topo)
    print(topoMap.value(-122,45))
    topoMap.profile(-131+360,43.0,-118+360,43.0)

    # age = np.random.randint(10,60,12)
    # Name = [f'Student{i:02d}' for i in range(1,13)]
    # grade = np.random.rand(12)*60

    # savetxt('test.txt',age,Name,grade)
    # data = loadtxt('test.txt')
    # savetxt2('test.txt',{'Name':Name,'age':age,'grade':grade},fmt=['<10','<4','<8.2f'],
    #          header=True,header_fmt=['<10','<4','<8'])
    # data = loadtxt2('test.txt')
    pass