import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from geographiclib.geodesic import Geodesic
from matplotlib.colors import LinearSegmentedColormap

labelFont = {'weight':'normal','size':14}
titleFont = {'weight':'bold','size':16}

import os,pycpt;
cvcpt = pycpt.load.gmtColormap(os.path.dirname(__file__)+'/cv_original.cpt')
rbcpt = pycpt.load.gmtColormap(os.path.dirname(__file__)+'/rainbow.cpt').reversed()

coastlineData = '/home/ayu/SubdInv/models/ne_10m_coastline/ne_10m_coastline'
physioData = os.path.dirname(__file__)+'/physio/physio'

def plotLocalBase(minlon, maxlon, minlat, maxlat,resolution='c',coastline=None,figwidth=None,ax=None,
             dlat=5.0, dlon=5.0, topo=None):
    ''' Plot base map with country, state boundaries '''
    rsphere = (6378137.00,6356752.3142)
    lon_centre, lat_centre = (minlon+maxlon)/2, (minlat+maxlat)/2
    distEW = Geodesic.WGS84.Inverse(minlat, minlon, minlat, maxlon)['s12']
    distNS = Geodesic.WGS84.Inverse(minlat, minlon, maxlat, minlon)['s12']
    
    if ax is None:
        figwidth = 6.0 if figwidth is None else figwidth
    else:
        figwidth = None

    if figwidth is not None:
        fig = plt.figure(figsize=[figwidth,figwidth/distEW*distNS])
    elif ax is not None:
        plt.sca(ax)
        fig = plt.gcf()
    else:
        print('Maybe something goes wrong in plotBase, please check!')
        fig = plt.figure()

    m = Basemap(width=distEW, height=distNS, rsphere=rsphere, resolution=resolution, projection='lcc', 
        lat_1=minlat, lat_2=maxlat, lon_0=lon_centre, lat_0=lat_centre)
    
    if topo is None:
        if coastline is not None:
            m.readshapefile(coastline,'coastline',linewidth=0.5)
        else:
            m.drawcoastlines()
        m.drawcountries(linewidth=1)
        m.drawparallels(np.arange(minlat,maxlat,dlat), labels=[1,0,0,0])
        m.drawmeridians(np.arange(minlon,maxlon,dlon), labels=[0,0,0,1])
        m.drawstates(color='k', linewidth=0.5,linestyle='solid')
        m.readshapefile(physioData,'physio',linewidth=0.25)
        return (fig,m)
    elif topo is True:
        m.etopo(scale=1.0)
        m.drawparallels(np.arange(minlat,maxlat,dlat), labels=[1,0,0,0])
        m.drawmeridians(np.arange(minlon,maxlon,dlon), labels=[0,0,0,1])
        return (fig,m)
    else:
        raise ValueError('Plot with topo data has not been done yet!')

def plotGlobalBase(figwidth=None,ax=None,resolution='l'):
    if ax is None:
        figwidth = 12.0 if figwidth is None else figwidth
    else:
        figwidth = None
    if figwidth is not None:
        fig = plt.figure(figsize=[figwidth,figwidth/2])
    elif ax is not None:
        plt.sca(ax)
        fig = plt.gcf()
    else:
        print('Maybe something goes wrong in plotBase, please check!')
        fig = plt.figure()
    m = Basemap(resolution='c')
    m.drawcoastlines(linewidth=1)
    return fig,m

def plotGlobalCart(type='Plate',**kwargs):
    crsPlate = ccrs.PlateCarree()
    if type == 'Plate':
        crs = ccrs.PlateCarree()
        plt.figure(figsize=[15,6])
    elif type == 'Azimuth':
        lat,lon = kwargs['central_latitude'],kwargs['central_longitude']
        crs = ccrs.AzimuthalEquidistant(central_longitude=lon, central_latitude=lat)
        plt.figure()
    ax = plt.axes(projection=crs)
    ax.set_global()
    ax.coastlines()
    return ax,crsPlate

def colorCycle():
    fig = plt.figure()
    ax = fig.add_subplot()
    cycle = ax._get_lines.prop_cycler
    cl = [next(cycle) for _ in range(10)]
    cl = [c['color'] for c in cl]
    plt.close()
    return cl

def cpt2cmap(cptfile,name='NewColorMap',N=256):
    ''' Generate linearly changed colormap from cpt file, B,F,N lines are ignored 
    Attention: when converted into cmap, the color palette table will be resampled. It is ok
    when cpt is continuous, however there may be a color jump sometime, like sea level when plot
    topo graph. The cmap struct will have N color ranges with N+1 boundaries. The color is
    resampled at np.linspace(0,1,N) and boundaries are np.linspace(0,1,N+1), after normalizing 
    z-value in cpt to [0,1]. Each color range is half-closed interval [a,b) except last one, which
    is closed interval [a,b].
    So if you want to get the color jump correctly, you can choose appropriate N here to make sure
    make sure one of the range boundaries locates at normalized jump point of cpt file. And then
    choose vmin and vmax when plotting to make sure the value for plot is also at that point. For example:
    if cpt is for topo plot, min = -5000, max = 8000, jump at z = 0
    N = 260, vmax/(-vmin) = 8/5 would be a good set of parameters.
    '''
    cptmat = []
    with open(cptfile) as f:
        for line in f:
            line = line.strip()
            if (not line) or (line[0] in ('#','B','F','N')):
                continue
            else:
                cptmat.append([float(i) for i in line.split()])

    zcheckb = np.array([cptmat[i][0] for i in range(len(cptmat))])
    zchecke = np.array([cptmat[i][4] for i in range(len(cptmat))])
    if np.any(zchecke - zcheckb < 0) or np.any(zchecke[:-1] - zcheckb[1:] != 0):
        raise ValueError('Illegal .cpt file! The z-values should increase monotonically.')

    red     = []
    blue    = []
    green   = []
    alpha   = []
    for i in range(len(cptmat)):
        zvalue = (cptmat[i][0]-cptmat[0][0])/(cptmat[-1][0]-cptmat[0][0])
        if i == 0:
            red.append(     (zvalue,0,cptmat[i][1]/255.) )
            green.append(   (zvalue,0,cptmat[i][2]/255.) )
            blue.append(    (zvalue,0,cptmat[i][3]/255.) )
            alpha.append(   (zvalue,1.0,1.0) )
        else:
            red.append(     (zvalue,cptmat[i-1][5]/255.,cptmat[i][1]/255.) )
            green.append(   (zvalue,cptmat[i-1][6]/255.,cptmat[i][2]/255.) )
            blue.append(    (zvalue,cptmat[i-1][7]/255.,cptmat[i][3]/255.) )
            alpha.append(   (zvalue,1.0,1.0) )
    red.append(     (zvalue,cptmat[i][5]/255.,1.0) )
    green.append(   (zvalue,cptmat[i][6]/255.,1.0) )
    blue.append(    (zvalue,cptmat[i][7]/255.,1.0) )
    alpha.append(   (zvalue,1.0,1.0) )
    cdict = {'red':red, 'green':green, 'blue':blue, 'alpha':alpha}
    return LinearSegmentedColormap(name,cdict,N)

def mkcmap(inM,name='NewColorMap',N=256):
    ''' 
    Generate linearly changed colormap, the larger N, the smoother the change would be.
    inM is a list of [z-value, R, G, B, (alpha)]
    '''
    inM = np.array(inM)
    n = inM.shape[0]
    if inM.shape[1] == 4:
        inM = np.append(inM,np.ones([n,1]),axis=1)
    inM[:,0] = (inM[:,0] - inM[0,0]) / (inM[-1,0] - inM[0,0])
    red     = []
    blue    = []
    green   = []
    alpha   = []
    for i in range(n):
        red.append((    inM[i,0],inM[i,1],inM[i,1]))
        green.append((  inM[i,0],inM[i,2],inM[i,2]))
        blue.append((   inM[i,0],inM[i,3],inM[i,3]))
        alpha.append((  inM[i,0],inM[i,4],inM[i,4]))
    cdict = {'red':red, 'green':green, 'blue':blue, 'alpha':alpha}
    return LinearSegmentedColormap(name,cdict,N)


# fig, axes = plt.subplots(2,3,figsize=[12,8])
# plt.subplots_adjust(wspace=0.25,hspace=0.3,left=0.08,right=0.92,bottom=0.08)



if __name__ == '__main__':
    import matplotlib.pyplot as plt
    plt.ion()
    print('Testing ...')
    plotLocalBase(-172.0, -122.0, 52.0, 72.0,resolution='l')
    fig,m = plotLocalBase(-130,-110,30,50,resolution='l')
    # fig,m = plotLocalBase(-130,-110,30,50,coastline=coastlineData,topo='True')
    # plotGlobalBase()
