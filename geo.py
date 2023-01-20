import numpy as np
from geographiclib.geodesic import Geodesic
from shapely.geometry import Polygon,Point

# assuming ellipsoid
def region_near_path(lat1,lon1,lat2,lon2,d,N=100,polyType='shapely',half=False): # d in meter
    l = Geodesic.WGS84.InverseLine(lat1,lon1,lat2,lon2)
    poly = []
    for v in np.arange(0,d,20):
        p1 = Geodesic.WGS84.Direct(lat1,lon1, l.azi1-90, v)
        p2 = Geodesic.WGS84.Direct(lat1,lon1, l.azi1+90, v)
        poly.append((p1['lon2'],p1['lat2']))
        poly.insert(0,(p2['lon2'],p2['lat2']))
    for p0 in [l.Position(v) for v in np.linspace(0,l.s13,N)]:
        p1 = Geodesic.WGS84.Direct(p0['lat2'],p0['lon2'], p0['azi2']-90, d)
        p2 = Geodesic.WGS84.Direct(p0['lat2'],p0['lon2'], p0['azi2']+90, d)
        poly.append((p1['lon2'],p1['lat2']))
        poly.insert(0,(p2['lon2'],p2['lat2']))
    for v in np.arange(0,d,20)[::-1]:
        p1 = Geodesic.WGS84.Direct(p0['lat2'],p0['lon2'], p0['azi2']-90, v)
        p2 = Geodesic.WGS84.Direct(p0['lat2'],p0['lon2'], p0['azi2']+90, v)
        poly.append((p1['lon2'],p1['lat2']))
        poly.insert(0,(p2['lon2'],p2['lat2']))
    if polyType == 'shapely':
        return Polygon(poly)
    elif polyType == 'matplotlib':
        import matplotlib as mpl
        return mpl.patches.Polygon(poly)

def region_near_point(lat0,lon0,d,N=300): # d in meter
    poly = []
    for p1 in [Geodesic.WGS84.Direct(lat0,lon0,azi,d) for azi in np.linspace(0,360)]:
        poly.append((p1['lon2'],p1['lat2']))
    import matplotlib as mpl
    return mpl.patches.Polygon(poly)

def region_along_path(path,d):
    def _interp(path,interp=1):
        pathInterp = []
        for i in range(len(path)-1):
            lon1,lat1 = path[i]; lon2,lat2 = path[i+1]
            l = Geodesic.WGS84.InverseLine(lat1,lon1,lat2,lon2)
            for d in np.linspace(0,l.s13/1e3,int(l.s13/1e3//interp)+2)[:-1]:
                p = Geodesic.WGS84.Direct(lat1,lon1, l.azi1, d*1e3)
                pathInterp.append((p['lon2'],p['lat2']))
        pathInterp.append(path[-1])
        pathInterp = np.asarray(pathInterp)
        return pathInterp
    def _region_along_line(p1,p2,d,interp=1e3,ending=None):
        lon1,lat1,lon2,lat2 = p1[0],p1[1],p2[0],p2[1]
        l = Geodesic.WGS84.InverseLine(lat1,lon1,lat2,lon2)
        poly = []
        for v in np.arange(0,d,d/3):
            p1 = Geodesic.WGS84.Direct(lat1,lon1, l.azi1-90, v)
            p2 = Geodesic.WGS84.Direct(lat1,lon1, l.azi1+90, v)
            poly.append((p1['lon2'],p1['lat2']))
            poly.insert(0,(p2['lon2'],p2['lat2']))
        for p0 in [l.Position(v) for v in np.linspace(0,l.s13,int(l.s13//interp)+2)]:
            p1 = Geodesic.WGS84.Direct(p0['lat2'],p0['lon2'], p0['azi2']-90, d)
            p2 = Geodesic.WGS84.Direct(p0['lat2'],p0['lon2'], p0['azi2']+90, d)
            poly.append((p1['lon2'],p1['lat2']))
            poly.insert(0,(p2['lon2'],p2['lat2']))
        for v in np.arange(0,d,d/3)[::-1]:
            p1 = Geodesic.WGS84.Direct(p0['lat2'],p0['lon2'], p0['azi2']-90, v)
            p2 = Geodesic.WGS84.Direct(p0['lat2'],p0['lon2'], p0['azi2']+90, v)
            poly.append((p1['lon2'],p1['lat2']))
            poly.insert(0,(p2['lon2'],p2['lat2']))
        return poly
    
    from shapely.geometry import Polygon,LineString
    from shapely.ops import unary_union
    # path=_interp(path)
    polys = [Polygon(_region_along_line(p1,p2,d)) for p1,p2 in zip(path[:-1],path[1:])]
    poly = unary_union(polys)
    if type(poly.boundary) == LineString:
        return np.asarray(poly.boundary.xy).T
    else:
        return np.asarray(poly.exterior.xy).T


# assuming sphere
def wgs2vec(lon,lat):
    return np.array([
        np.cos(np.radians(lat))*np.cos(np.radians(lon)),
        np.cos(np.radians(lat))*np.sin(np.radians(lon)),
        np.sin(np.radians(lat))
        ])
def vec2wgs(vec):
    eps = 1e-6
    x,y,z = vec/np.sqrt((vec**2).sum())
    if abs(1-abs(z))<eps:
        return 0,90*np.sign(z)
    lat = np.rad2deg(np.arcsin(z))
    lon = np.rad2deg(np.arctan2(y,x))
    return lon,lat
def point_project_to_profile(lon1,lat1,lon2,lat2,lon0,lat0,R=6371):
    vec1,vec2,vec0 = wgs2vec(lon1,lat1),wgs2vec(lon2,lat2),wgs2vec(lon0,lat0)
    vecV = np.cross(vec1,vec2)
    vecP = vec0-np.inner(vecV,vec0)/np.inner(vecV,vecV)*vecV
    lon,lat = vec2wgs(vecP)
    dist = np.arccos(np.inner(vec0,vecP)/np.sqrt(np.inner(vec0,vec0)*np.inner(vecP,vecP)))*R
    return lon,lat,dist
def sphereDist(lon1,lat1,lon2,lat2,R=6371):
    return np.arccos(np.inner(wgs2vec(lon1,lat1),wgs2vec(lon2,lat2)))*R
    

from Triforce.utils import GeoGrid,GeoMap

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


if __name__ == '__main__':
    from Triforce.pltHead import *
    
    ''' region_near_path '''
    lat1,lon1,lat2,lon2 = 44.57461, -130.46451, 43.07297, -126.62971
    l = Geodesic.WGS84.InverseLine(lat1,lon1,lat2,lon2)
    poly1 = region_near_path(lat1,lon1,lat2,lon2,80e3,N=100)
    # poly2 = region_near_point(lat1,lon1,d=30e3)
    plt.figure()
    plt.plot((lon1,lon2),(lat1,lat2))
    x,y = poly1.boundary.xy
    plt.plot(x,y,'-o')
    # x,y = poly2.boundary.xy
    # plt.plot(x,y,'--')
    # poly1.contains(Point(-129,42.5))

    ''' region_along_path '''
    import shapefile
    from Triforce.customPlot import plotLocalBase
    PA_JF_Division = shapefile.Reader("/home/ayu/Projects/Cascadia/Models/Plates/PB2002_boundaries.shp").shapes()[66]
    poly1 = region_along_path([(-129.572,46.4015),(-128.94,48.3236)],d=80e3)
    poly2 = region_along_path(PA_JF_Division.points[4:20],d=80e3)
    fig,m = plotLocalBase(-132,-121,39.5,50,resolution='i',
                        gridlines={'dlat':3,'dlon':4,'lat0':42,'lon0':-130},
                        coastline=False,countries=False,states=False,plateBoundary=True)
    x,y = poly1.T
    m.plot(x,y,latlon=True)
    x,y = poly2.T
    m.plot(x,y,latlon=True)
    plt.figure()
    x,y = np.asarray(PA_JF_Division.points).T;plt.plot(x,y)
    x,y = poly1.T;plt.plot(x,y)
    x,y = poly2.T;plt.plot(x,y)