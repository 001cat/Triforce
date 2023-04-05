import numpy as np
from geographiclib.geodesic import Geodesic

class GeomMaker():
    '''All units are in m'''
    def __init__(self,R=6371) -> None:
        self.wgs = Geodesic.WGS84

    def _brush_dot(self,lon0,lat0,d,N=300):
        poly = []
        for p1 in [Geodesic.WGS84.Direct(lat0,lon0,azi,d) for azi in np.linspace(0,360)]:
            poly.append((p1['lon2'],p1['lat2']))
        return poly
    
    def _brush_line(self,p1,p2,d,interp=1e3):
        lon1,lat1,lon2,lat2 = p1[0],p1[1],p2[0],p2[1]
        l = self.wgs.InverseLine(lat1,lon1,lat2,lon2)
        poly = []
        for v in np.arange(0,d,d/3):
            p1 = self.wgs.Direct(lat1,lon1, l.azi1-90, v)
            p2 = self.wgs.Direct(lat1,lon1, l.azi1+90, v)
            poly.append((p1['lon2'],p1['lat2']))
            poly.insert(0,(p2['lon2'],p2['lat2']))
        for p0 in [l.Position(v) for v in np.linspace(0,l.s13,int(l.s13//interp)+2)]:
            p1 = self.wgs.Direct(p0['lat2'],p0['lon2'], p0['azi2']-90, d)
            p2 = self.wgs.Direct(p0['lat2'],p0['lon2'], p0['azi2']+90, d)
            poly.append((p1['lon2'],p1['lat2']))
            poly.insert(0,(p2['lon2'],p2['lat2']))
        for v in np.arange(0,d,d/3)[::-1]:
            p1 = self.wgs.Direct(p0['lat2'],p0['lon2'], p0['azi2']-90, v)
            p2 = self.wgs.Direct(p0['lat2'],p0['lon2'], p0['azi2']+90, v)
            poly.append((p1['lon2'],p1['lat2']))
            poly.insert(0,(p2['lon2'],p2['lat2']))
        return poly

    def _brush_path(self,path,d,interp=1e3):
        '''Does not work well, test only'''
        from shapely.geometry import Polygon,LineString
        from shapely.ops import unary_union
        polys = [Polygon(self._brush_line(p1,p2,d,interp)) for p1,p2 in zip(path[:-1],path[1:])]
        poly = unary_union(polys)
        if type(poly.boundary) == LineString:
            return np.asarray(poly.boundary.xy).T
        else:
            return np.asarray(poly.exterior.xy).T

class GeomMaker_Sphere():
    '''DEBUGING, DO NOT USE NOW!!!!'''
    def __init__(self,R=6371) -> None:
        self.R = R

    @staticmethod
    def wgs2vec(lon,lat):
        return np.array([
            np.cos(np.radians(lat))*np.cos(np.radians(lon)),
            np.cos(np.radians(lat))*np.sin(np.radians(lon)),
            np.sin(np.radians(lat))
            ])
    @staticmethod
    def vec2wgs(vec):
        eps = 1e-6
        x,y,z = vec/np.sqrt((vec**2).sum())
        if abs(1-abs(z))<eps:
            return 0,90*np.sign(z)
        lat = np.rad2deg(np.arcsin(z))
        lon = np.rad2deg(np.arctan2(y,x))
        return lon,lat
    @staticmethod
    def decomponent(vec,vec0):
        vecP = np.inner(vec,vec0)/np.inner(vec0,vec0)*vec0 # parallel to vec0
        vecV = vec-vecP # perpendicular (vertical) to vec0
        return vecP,vecV
    @staticmethod
    def norm2(vec):
        return vec/np.sqrt(np.inner(vec,vec))

    def dist(self,lon1,lat1,lon2,lat2):
        return np.arccos(np.inner(self.wgs2vec(lon1,lat1),self.wgs2vec(lon2,lat2)))*self.R
    
    def azi(self,lon1,lat1,lon2,lat2):
        if abs(lat1) == 90 or abs(lat2) == 90 or (lat1+lat2==0 and abs(lon1-lon2)==180):
            raise ValueError('Antipoles are found!')
        if (lat1 == lat2 and lon1 == lon2):
            raise ValueError('Points overlaped')
        
        vec0,vec1,vec2 = np.array([0,0,1]),self.wgs2vec(lon1,lat1),self.wgs2vec(lon2,lat2)

        _,vecV10 = self.decomponent(vec1,vec0); vecV10/=np.sqrt(np.inner(vecV10,vecV10))
        _,vecV12 = self.decomponent(vec1,vec2); vecV12/=np.sqrt(np.inner(vecV12,vecV12))
        azi12 = np.arccos(np.inner(vecV10,vecV12))/np.pi*180
        if np.inner(np.cross(vecV10,vecV12),vec1) > 0:
            azi12 = -azi12
        
        _,vecV20 = self.decomponent(vec2,vec0); vecV20/=np.sqrt(np.inner(vecV20,vecV20))
        _,vecV21 = self.decomponent(vec2,vec1); vecV21/=np.sqrt(np.inner(vecV21,vecV21))
        azi21 = np.arccos(np.inner(vecV20,vecV21))/np.pi*180
        # from IPython import embed; embed()
        if np.inner(np.cross(vecV20,vecV21),vec2) > 0:
            azi21 = -azi21
        azi21 -= 180*np.sign(azi21) # to keep consistent with geographiclib.geodesics
        
        return azi12,azi21

    def forward(self,lon1,lat1,azi,dist):
        from scipy.spatial.transform import Rotation as R
        def rotate_clockwise(vec0,vec1,theta):
            return R.from_rotvec(-theta*self.norm2(vec0),degrees=True).apply(vec1)
        vec0 = np.array([0,0,1]); vec1 = self.wgs2vec(lon1,lat1)
        vecV = np.cross(vec0,vec1)
        vecV = rotate_clockwise(vec1,vecV,azi)
        vec2 = rotate_clockwise(vecV,vec1,dist/self.R/np.pi*180)
        lon2,lat2 = self.vec2wgs(vec2)
        return lon2,lat2

    def point_project_to_profile(self,lon1,lat1,lon2,lat2,lon0,lat0):
        vec1,vec2,vec0 = self.wgs2vec(lon1,lat1),self.wgs2vec(lon2,lat2),self.wgs2vec(lon0,lat0)
        vecV = np.cross(vec1,vec2)
        vecP = vec0-np.inner(vecV,vec0)/np.inner(vecV,vecV)*vecV
        lon,lat = self.vec2wgs(vecP)
        dist = self.dist(lon,lat,lon0,lat0)
        return lon,lat,dist

    def _brush_line(self,p1,p2,d):
        from shapely.geometry import Polygon
        lon1,lat1,lon2,lat2 = p1[0],p1[1],p2[0],p2[1]
        azi12,azi21 = self.azi(lon1,lat1,lon2,lat2)
        p1 = self.forward(lon1,lat1,azi12-90,d)
        p2 = self.forward(lon1,lat1,azi12+90,d)
        p3 = self.forward(lon2,lat2,azi21+90,d)
        p4 = self.forward(lon2,lat2,azi21-90,d)
        return [p1,p2,p3,p4]
        



if __name__ == '__main__':
    from Triforce.pltHead import *
    from Triforce.cartopy import plotLocalCart,crsPlate
    

    geomMaker = GeomMaker_Sphere()

    print(geomMaker.azi(0,80,40,80))
    print(Geodesic.WGS84.Inverse(80,0,80,40))

    from scipy.spatial.transform import Rotation as R
    r = R.from_rotvec(30*np.array([0,0,1]),degrees=True)
    r.apply([0,1,0])



    polygon1 = GeomMaker_Sphere()._brush_line((-129.453,46.487),(-128.906,47.709),100)
    polygon2 = GeomMaker()._brush_line((-129.453,46.487),(-128.906,47.709),100e3)

    minlon,maxlon,minlat,maxlat = -140,-110,35,55
    ax,_ = plotLocalCart(minlon,maxlon,minlat,maxlat)
    lons = [p[0] for p in polygon1]
    lats = [p[1] for p in polygon1]
    ax.plot(lons,lats,transform=crsPlate)
    lons = [p[0] for p in polygon2]
    lats = [p[1] for p in polygon2]
    ax.plot(lons,lats,transform=crsPlate)



    path = [(-127.701,40.571),(-127.557,41.560),(-127.217,41.708),(-126.626,43.023)]
    polygon3 = GeomMaker()._brush_path(path,100e3)
    lons = [p[0] for p in polygon3]
    lats = [p[1] for p in polygon3]
    ax.plot(lons,lats,transform=crsPlate)