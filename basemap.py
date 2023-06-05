import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from typing import Tuple
from mpl_toolkits.basemap import Basemap
from geographiclib.geodesic import Geodesic

# coastlineData = '/home/ayu/Packages/ne_10m_coastline/ne_10m_coastline'
physioData = os.path.dirname(__file__)+'/physio/physio' #https://water.usgs.gov/GIS/dsdl/physio_shp.zip


''' Plotting Geographic Maps '''
def plotLocalBase(minlon, maxlon, minlat, maxlat, projection='merc', resolution='c',
                  ax=None, topo=False, coastline=False, countries=True, states=True, plateBoundary=True,
                  gridlines = {
                    'dlat':5.0,'lat0':None,'latLocation':[1,0,0,0],
                    'dlon':5.0,'lon0':None,'lonLocation':[0,0,0,1]
                    }) -> Tuple[mpl.figure.Figure,Basemap]:
    ''' Plot base map with country, state boundaries '''
    minlon,maxlon = minlon-360*(minlon>180),maxlon-360*(maxlon>180)
    rsphere = (6378137.00,6356752.3142)
    lon_centre, lat_centre = (minlon+maxlon)/2, (minlat+maxlat)/2
    distEW = Geodesic.WGS84.Inverse(minlat, minlon, minlat, maxlon)['s12']
    distNS = Geodesic.WGS84.Inverse(minlat, minlon, maxlat, minlon)['s12']


    if gridlines:
        gridlines_tmp = gridlines
        gridlines = {
            'dlat':5.0,'lat0':None,'latLocation':[1,0,0,0],
            'dlon':5.0,'lon0':None,'lonLocation':[0,0,0,1]
        }
        gridlines.update(gridlines_tmp)

    if ax is None:
        figwidth = 6.0
        fig = plt.figure(figsize=[figwidth,figwidth/distEW*distNS])
        ax = plt.gca()
    else:
        plt.sca(ax); fig = plt.gcf()

    if projection == 'lcc':
        m = Basemap(width=distEW, height=distNS, rsphere=rsphere, resolution=resolution, projection='lcc', 
            lat_1=minlat, lat_2=maxlat, lon_0=lon_centre, lat_0=lat_centre)
    elif projection in ['merc','mill']:
        m = Basemap(projection=projection, llcrnrlat=minlat, urcrnrlat=maxlat, 
                    llcrnrlon=minlon, urcrnrlon=maxlon,
                    lat_ts=(minlat+maxlat)/2,resolution=resolution)
    else:
        raise ValueError('Not supported yet.')

    ax.set_facecolor([0.9]*3)
    [ax.spines[l].set_linewidth(0.25) for l in ['top','bottom','left','right']]
    m.fillcontinents()

    if plateBoundary:
        m.readshapefile(os.path.dirname(__file__)+"/PB2002/PB2002_boundaries",
                    'PB2002_boundaries',linewidth=2.0,color='orange')

    if topo is True:
        m.etopo(scale=1.0)

    if coastline == False:
        pass
    elif coastline == True:
        m.drawcoastlines()
    elif coastline == 'accurate':
        m.readshapefile(coastline,'coastline',linewidth=0.5)
    else:
        raise ValueError('Customized coastline data has not been supported yet!')

    if countries:
        m.drawcountries(linewidth=1)

    if states:
        m.drawstates(color='k', linewidth=0.5,linestyle='solid')
    
    if gridlines:
        lat0_tick = gridlines['lat0'] or minlat
        lon0_tick = gridlines['lon0'] or minlon
        dlat,dlon = gridlines['dlat'],gridlines['dlon']
        PM = {}
        PM.update(m.drawparallels(np.arange(lat0_tick,maxlat+dlat/2,dlat),labels=gridlines['latLocation']))
        PM.update(m.drawmeridians(np.arange(lon0_tick,maxlon+dlon/2,dlon),labels=gridlines['lonLocation']))
        for v in PM.values():
            v[0][0].set_alpha(0.2)
    return (fig,m)

def plotGlobalBase(figwidth=None,ax=None,resolution='l') -> Tuple[mpl.figure.Figure,Basemap]:
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
    # m.drawcoastlines(linewidth=1)
    m.fillcontinents()
    return fig,m

def plotBasemap_Cascadia(loc=(-132,-121,39.5,50,4,3,-130,42),ax=None,paLable=[1,0,0,0],meLable=[0,0,0,1]):
    minlon,maxlon,minlat,maxlat,dlon,dlat,lon0_tick,lat0_tick = loc
    gridlines = {
        'dlat':dlat,'lat0':lat0_tick,'latLocation':paLable,
        'dlon':dlon,'lon0':lon0_tick,'lonLocation':meLable
    }
    fig,m = plotLocalBase(minlon,maxlon,minlat,maxlat,resolution='i',ax=ax,
                          gridlines=gridlines,coastline=False,countries=False,states=False,
                          plateBoundary=True)
    return fig,m
