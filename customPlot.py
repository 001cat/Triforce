import os
import numpy as np
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from geographiclib.geodesic import Geodesic
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

labelFont = {'weight':'normal','size':14}
titleFont = {'weight':'bold','size':16}

# import os,pycpt;
# cvcpt = pycpt.load.gmtColormap(os.path.dirname(__file__)+'/cv_original.cpt')
# rbcpt = pycpt.load.gmtColormap(os.path.dirname(__file__)+'/rainbow.cpt').reversed()

coastlineData = '/home/ayu/SubdInv/models/ne_10m_coastline/ne_10m_coastline'
physioData = os.path.dirname(__file__)+'/physio/physio' #https://water.usgs.gov/GIS/dsdl/physio_shp.zip

def plotLocalBase(minlon, maxlon, minlat, maxlat,resolution='c',coastline=None,figwidth=None,ax=None,
             dlat=5.0, dlon=5.0, topo=None, projection='merc'):
    ''' Plot base map with country, state boundaries '''
    minlon,maxlon = minlon-360*(minlon>180),maxlon-360*(maxlon>180)
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

    if projection == 'lcc':
        m = Basemap(width=distEW, height=distNS, rsphere=rsphere, resolution=resolution, projection='lcc', 
            lat_1=minlat, lat_2=maxlat, lon_0=lon_centre, lat_0=lat_centre)
    elif projection in ['merc','mill']:
        m = Basemap(projection=projection, llcrnrlat=minlat, urcrnrlat=maxlat, 
                    llcrnrlon=minlon, urcrnrlon=maxlon,
                    lat_ts=(minlat+maxlat)/2,resolution=resolution)
    else:
        raise ValueError('Not supported yet.')
    m.readshapefile('/home/ayu/Projects/Cascadia/Models/Plates/PB2002_boundaries',
                'PB2002_boundaries',linewidth=2.0,color='orange')
    
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

def plotLocalCart(minlon, maxlon, minlat, maxlat, dlat=5.0, dlon=5.0, fig=None):
    crsPlate = ccrs.PlateCarree()
    # crs = ccrs.PlateCarree()
    crs = ccrs.Mercator()
    if fig is None:
        plt.figure()
    else:
        plt.figure(fig.number)
    ax = plt.axes(projection=crs)
    ax.set_extent((minlon, maxlon, minlat, maxlat))
    ax.coastlines()

    import matplotlib.ticker as mticker
    gl = ax.gridlines(crs=crsPlate, draw_labels=True,linewidth=1, color='gray', 
                      alpha=0.5, linestyle='--')
    gl.xlocator = mticker.FixedLocator(np.arange(minlon-360*(minlon>180),maxlon-360*(maxlon>180),dlon))
    gl.ylocator = mticker.FixedLocator(np.arange(minlat,maxlat,dlat))
    gl.xlabels_top      = False
    gl.xlabels_bottom   = True
    gl.ylabels_left     = True
    gl.ylabels_right    = False

    import cartopy.feature as cfeature
    states_provinces = cfeature.NaturalEarthFeature(
            category='cultural',  name='admin_1_states_provinces_lines',
            scale='50m', facecolor='none')
    ax.add_feature(states_provinces, edgecolor='black', zorder=10) 

    import cartopy.io.shapereader as shpreader
    import cartopy.feature as cfeature
    # Read shape file
    reader = shpreader.Reader("/home/ayu/Projects/Cascadia/Models/Plates/PB2002_boundaries.shp")
    # Filter for a specific country
    plateBoundary = [plate.geometry for plate in reader.records()]
    shape_feature = cfeature.ShapelyFeature(plateBoundary, crsPlate, facecolor=(1, 1, 1, 0), 
                                            edgecolor='orange', lw=2)
    ax.add_feature(shape_feature)


    # ax.set_xticks([230-360,240-360], crs=crsPlate)
    # ax.set_xticklabels([230-360,240-360], color='red', weight='bold')
    # ax.set_yticks([40,45,50], crs=crsPlate)
    # ax.set_yticklabels([40,45,50])
    # ax.yaxis.tick_right()

    return ax,crsPlate

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
    def interpColorList(listedCmap,N):
        R = np.linspace(listedCmap[0][0],listedCmap[1][0],N)
        G = np.linspace(listedCmap[0][1],listedCmap[1][1],N)
        B = np.linspace(listedCmap[0][2],listedCmap[1][2],N)
        A = np.linspace(listedCmap[0][3],listedCmap[1][3],N)
        return np.array([R,G,B,A]).T

    cptmat = []
    with open(cptfile) as f:
        for line in f:
            line = line.strip()
            if (not line) or (line[0] in ('#','B','F','N')):
                continue
            else:
                cptmat.append([float(i) for i in line.split()])
    cptmat = np.array(cptmat)
    zLs = cptmat[:,4] - cptmat[:,0]
    Ns = np.round(zLs/zLs.sum()*N)
    Ns[-1] += N-Ns.sum()
    colorList = np.zeros((N,4))
    normX = np.zeros(len(zLs)+1)
    normY = np.zeros(len(zLs)+1)
    count = 0
    for i in range(len(zLs)):
        colorList[count:count+int(Ns[i])] = interpColorList(
                                    [[cptmat[i,1]/255,cptmat[i,2]/255,cptmat[i,3]/255,1],
                                     [cptmat[i,5]/255,cptmat[i,6]/255,cptmat[i,7]/255,1]],
                                    N=int(Ns[i]))
        normX[i] = cptmat[i,0]
        normY[i] = Ns[:i].sum()/N
        count += int(Ns[i])
    normX[-1] = cptmat[-1,4]
    normY[-1] = 1
    cmap = ListedColormap(colorList)
    norm = CustomNorm(normX,normY)
    return cmap,norm

class CustomNorm(mpl.colors.Normalize):
    def __init__(self,x,y,clip=False):
        super().__init__(x[0],x[-1],clip)
        self.x,self.y = x,y
    def __call__(self,value,clip=None):
        if type(value) == np.ma.masked_array:
            return np.ma.masked_array(np.interp(value,self.x,self.y),mask=value.mask)
        else:
            return np.ma.masked_array(np.interp(value,self.x,self.y))

def addAxes(loc,hasTitle=True):
    ''' loc = [x,y,w,h] '''
    if hasTitle:
        border = [0.125,0.11,0.1,0.12]
    else:
        border = [0.125,0.11,0.1,0.08]
    x,y,w,h = loc
    dx1,dy1,dx2,dy2 = border
    rect = [x+w*dx1, y+h*dy1, w*(1-dx1-dx2), h*(1-dy1-dy2)]
    return plt.axes(rect)

def addCAxes(ax,location='right',size=0.05,pad=0.12):
    bbox = ax.get_position()
    x,y,w,h = bbox.x0,bbox.y0,bbox.width,bbox.height
    if location == 'right':
        pad = 0.05 if pad is None else pad
        rect = [x+w+w*pad,y,w*size,h]
    elif location == 'bottom':
        pad = 0.12 if pad is None else pad
        rect = [x,y-w*pad-w*size,w,h*size]
    else:
        raise ValueError('Not ready yet')
    return plt.axes(rect)

def moveAxes(ax,move):
    bbox = ax.get_position()
    x,y,w,h = bbox.x0,bbox.y0,bbox.width,bbox.height
    if move[0] == 'R':
        x = x+float(move[1:])
    elif move[0] == 'L':
        x = x-float(move[1:])
    elif move[0] == 'U':
        y = y+float(move[1:])
    elif move[0] == 'D':
        y = y-float(move[1:])
    else:
        raise ValueError('Wrong move command, should start with U/D/R/L')
    ax.set_position([x,y,w,h])

# fig, axes = plt.subplots(2,3,figsize=[12,8])
# plt.subplots_adjust(wspace=0.25,hspace=0.3,left=0.08,right=0.92,bottom=0.08)

cvcpt = cpt2cmap(os.path.dirname(__file__)+'/cv_original.cpt','cv')[0]
rbcpt = cpt2cmap(os.path.dirname(__file__)+'/rainbow.cpt','rainbow')[0].reversed()


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    plt.ion()
    # print('Testing ...')
    # plotLocalBase(-172.0, -122.0, 52.0, 72.0,resolution='l')
    # fig,m = plotLocalBase(-130,-110,30,50,resolution='l')
    # # fig,m = plotLocalBase(-130,-110,30,50,coastline=coastlineData,topo=True)
    # fig,m = plotGlobalBase()
    # ax,crsPlate = plotGlobalCart()
    # ax,crsPlate = plotLocalCart(-130,-110,30,50)


    # plt.figure()
    # addAxes([0,0,1/2,1])
    # addAxes([1/2,0,1/2,1],hasTitle=False)

    # cmap,norm = cpt2cmap('GMT_globe.cpt')
    # from netCDF4 import Dataset
    # with Dataset('/home/ayu/Packages/etopo2/ETOPO2v2g_f4.nc') as etopo:
    # 	lon = etopo.variables['x'][:]
    # 	lat = etopo.variables['y'][:]
    # 	topo = etopo.variables['z'][:]
    # lat_bnd, lon_bnd = [52, 72], [-172, -122]
    # lat_ind = np.where((lat > lat_bnd[0]) & (lat < lat_bnd[1]))
    # lon_ind = np.where((lon > lon_bnd[0]) | (lon < lon_bnd[1]))
    # topo = topo[lat_ind[0],:]
    # topo = topo[:,lon_ind[0]]
    # lon  = lon[lon_ind]
    # lat  = lat[lat_ind]
    # XX,YY = np.meshgrid(lon,lat)
    # fig,m = plotLocalBase(-172.0, -122.0, 52.0, 72.0,resolution='l')
    # m.pcolormesh(XX, YY, topo, shading='gouraud',cmap=cmap,norm=norm,latlon=True)

    # from mpl_toolkits.axes_grid1 import make_axes_locatable
    # plt.figure()
    # ax = plt.subplot()
    # im = ax.imshow(np.arange(100).reshape((10, 10)))
    # # create an axes on the right side of ax. The width of cax will be 5%
    # # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes("right", size="5%", pad=0.50)
    # plt.colorbar(im, cax=cax)
    # # plt.colorbar(orientation='horizontal',fraction=0.1,aspect=50,pad=0.08)


    # plt.figure()
    # ax = addAxes([0,0,1,1])
    # moveAxes(ax,'U0.06')
    # im = ax.imshow(np.arange(100).reshape((10, 10)))
    # cax = addCAxes(ax,'bottom')
    # plt.colorbar(im,cax=cax,orientation="horizontal")
