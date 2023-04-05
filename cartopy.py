import os
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from typing import Tuple
from geographiclib.geodesic import Geodesic



import matplotlib.ticker as mticker
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from cartopy.mpl.geoaxes import GeoAxes
from cartopy.crs import PlateCarree

crsPlate = ccrs.PlateCarree()
crsMerc  = ccrs.Mercator()
def _createGeoAxes(projection) -> GeoAxes:
    # wrap function for using pylint
    return plt.axes(projection=projection)
def plotLocalCart(minlon, maxlon, minlat, maxlat, projection='merc', 
                  ax=None, topo=False, coastline=False, countries=True, states=True, plateBoundary=True, 
                  gridlines = {
                    'dlat':5.0,'lat0':None,'latLocation':[1,0,0,0],
                    'dlon':5.0,'lon0':None,'lonLocation':[0,0,0,1]
                    }) -> Tuple[GeoAxes,PlateCarree]:
    crsPlate = ccrs.PlateCarree()
    if projection == 'merc':
        crs = ccrs.Mercator()
    else:
        raise ValueError('Not supported yet.')

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
        ax = _createGeoAxes(projection=crs)
    else:
        plt.sca(ax); fig = plt.gcf()
    ax.set_extent((minlon, maxlon, minlat, maxlat))

    ax.set_facecolor([0.9]*3)
    [ax.spines[l].set_linewidth(0.25) for l in ['top','bottom','left','right']]
    # m.fillcontinents()
    ax.add_feature(cfeature.NaturalEarthFeature(
        category='physical', name='land', scale='50m', edgecolor='face', facecolor=[0.8]*3),zorder=0)

    if plateBoundary:
        reader = shpreader.Reader("/home/ayu/Projects/Cascadia/Models/Plates/PB2002_boundaries.shp")
        plateBoundary = [plate.geometry for plate in reader.records()]
        shape_feature = cfeature.ShapelyFeature(plateBoundary, crsPlate, facecolor=(1, 1, 1, 0), 
                                                edgecolor='orange', lw=2)
        ax.add_feature(shape_feature)

    if topo is True:
        etopo1 = os.path.dirname(__file__)+'/etopo1.jpg'
        if not os.path.exists(etopo1):
            raise ValueError(f'''Etopo background is not found! Please prepare it yourself 
            (https://github.com/matplotlib/basemap/blob/develop/packages/basemap_data/src/mpl_toolkits/basemap_data/etopo1.jpg)''')
        ax.imshow(plt.imread(os.path.dirname(__file__)+'/etopo1.jpg'),origin='upper',extent=(-180,180,-90,90),transform=crsPlate)
        # ax.stock_img()

    if coastline:
        ax.coastlines(resolution='50m')

    if countries:
        states_provinces = cfeature.NaturalEarthFeature(
            # category='cultural',  name='admin_0_boundary_lines_land',
            category='cultural',  name='admin_0_countries',
            scale='50m', facecolor='none')
        ax.add_feature(states_provinces, edgecolor='black', linewidth=1, zorder=10) 

    if states:
        states_provinces = cfeature.NaturalEarthFeature(
            category='cultural',  name='admin_1_states_provinces_lines',
            scale='50m', facecolor='none')
        ax.add_feature(states_provinces, edgecolor='black', linewidth=0.5, zorder=10) 


    if gridlines:
            lat0_tick = gridlines['lat0'] or minlat
            lon0_tick = gridlines['lon0'] or minlon
            dlat,dlon = gridlines['dlat'],gridlines['dlon']
            gl = ax.gridlines(crs=crsPlate, draw_labels=True,linewidth=1, color='k', alpha=0.2, linestyle=':')
            gl.xlocator     = mticker.FixedLocator(np.arange(lon0_tick,maxlon+dlon/2,dlon))
            gl.ylocator     = mticker.FixedLocator(np.arange(lat0_tick,maxlat+dlat/2,dlat))
            gl.top_labels      = gridlines['lonLocation'][2] == 1
            gl.bottom_labels   = gridlines['lonLocation'][3] == 1
            gl.left_labels     = gridlines['latLocation'][0] == 1
            gl.right_labels    = gridlines['latLocation'][1] == 1

    return ax,crsPlate

def plotGlobalCart(mapType='Plate',**kwargs):
    crsPlate = ccrs.PlateCarree()
    if mapType == 'Plate':
        crs = ccrs.PlateCarree()
        plt.figure(figsize=[15,6])
    elif mapType == 'Azimuth':
        lat,lon = kwargs['central_latitude'],kwargs['central_longitude']
        crs = ccrs.AzimuthalEquidistant(central_longitude=lon, central_latitude=lat)
        plt.figure()
    # ax = plt.axes(projection=crs)
    ax = _createGeoAxes(projection=crs)
    ax.set_global()
    ax.coastlines()
    return ax,crsPlate



# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# import matplotlib.pyplot as plt

# fig = plt.figure(figsize=(3, 3))

# ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertCylindrical())
# ax.set_extent([-20, 60, -40, 45])
# ax.add_feature(cfeature.LAND)
# ax.add_feature(cfeature.OCEAN)
# ax.add_feature(cfeature.COASTLINE)
# ax.add_feature(cfeature.BORDERS, linestyle=':')
# ax.add_feature(cfeature.LAKES, alpha=0.5)
# ax.add_feature(cfeature.RIVERS)
# ax.set_title('LambertCylindrical')

# plt.savefig('LambertCylindrical.png', dpi=100)