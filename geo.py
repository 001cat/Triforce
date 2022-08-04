import numpy as np
from geographiclib.geodesic import Geodesic
from shapely.geometry import Polygon,Point

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
    return Polygon(poly)

if __name__ == '__main__':
    from Triforce.pltHead import *
    
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