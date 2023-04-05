import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

labelFont = {'weight':'normal','size':14}
titleFont = {'weight':'bold','size':16}
text_bbox = {'fc': '0.8', 'pad': 2}

''' useful segments '''
# plt.errorbar(x,y,yerr,ls='None',marker='o',capsize=3,capthick=2,elinewidth=2)

# fig, axes = plt.subplots(2,3,figsize=[12,8])
# plt.subplots_adjust(wspace=0.25,hspace=0.3,left=0.08,right=0.92,bottom=0.08)
# text.set_path_effects([
#     path_effects.Stroke(linewidth=3,foreground='black'),
#     path_effects.Normal()
# ])

def labelSubplot(ax,label,x=0.03,y=0.97,kwargs={}):
    bbox = dict(boxstyle='round,pad=0.1', facecolor='white', alpha=1, linewidth=0.1)
    kwargsIn = dict(fontsize=18,fontweight='bold',ha='left',va='bottom',bbox=bbox)
    kwargsIn.update(kwargs)
    return ax.text(x,y,label,transform=ax.transAxes,**kwargsIn)


''' Axis modifications'''
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

def addCAxes(ax,location='right',size=0.05,pad=0.05):
    if type(ax) is list:
        axes = ax
        bounds = np.zeros((len(axes),4))
        for i,ax in enumerate(axes):
            bbox = ax.get_position()
            bounds[i,:2] = bbox.intervalx
            bounds[i,2:] = bbox.intervaly
        x,y,w,h = bounds[:,0].min(),bounds[:,2].min(),\
                  bounds[:,1].max()-bounds[:,0].min(),\
                  bounds[:,3].max()-bounds[:,2].min()
    else:
        bbox = ax.get_position()
        x,y,w,h = bbox.x0,bbox.y0,bbox.width,bbox.height
    if location == 'right':
        pad = 0.05 if pad is None else pad
        rect = [x+w+w*pad,y,w*size,h]
    elif location == 'bottom':
        pad = 0.05 if pad is None else pad
        rect = [x,y-h*pad-h*size,w,h*size]
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

def moveSpine(ax,move):
    bbox = ax.get_position()
    x0,x1 = bbox.intervalx
    y0,y1 = bbox.intervaly
    if move[0] == 'R':
        x1 += float(move[1:])
    elif move[0] == 'L':
        x0 += float(move[1:])
    elif move[0] == 'U':
        y1 += float(move[1:])
    elif move[0] == 'D':
        y0 += float(move[1:])
    else:
        raise ValueError('Wrong move command, should start with U/D/R/L')
    ax.set_position([x0,y0,x1-x0,y1-y0])

def scaleGeoAxes(ax,factor = 1):
    pos = ax.get_position()
    # ratio = pos.width/pos.height
    pos.x1 = (pos.x1-pos.x0)*factor+pos.x0
    pos.y1 = (pos.y1-pos.y0)*factor+pos.y0
    ax.set(position=pos)

''' colormap, cpt, norm '''
class CustomNorm(mpl.colors.Normalize):
    '''Please add cax.set_xscale('linear') to get a normal colorbar'''
    def __init__(self,x,y,clip=False):
        super().__init__(x[0],x[-1],clip)
        self.x,self.y = x,y
    def __call__(self,value,clip=None):
        if type(value) == np.ma.masked_array:
            return np.ma.masked_array(np.interp(value,self.x,self.y),mask=value.mask)
        else:
            return np.ma.masked_array(np.interp(value,self.x,self.y))
    def inverse(self,value):
        if np.iterable(value):
            if type(value) == np.ma.masked_array:
                return np.ma.masked_array(np.interp(value,self.y,self.x),mask=value.mask)
        return np.interp(value,self.y,self.x)

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
                cptmat.append([float(i) for i in line.replace('/',' ').split()])
    cptmat = np.array(cptmat)
    zLs = cptmat[:,4] - cptmat[:,0]

    # Ns = np.round(zLs/zLs.sum()*N)    # depracated, might result in negative Ns[-1]
    # Ns[-1] += N-Ns.sum()
    def calNs(zLs,N):
        Ns = np.zeros(zLs.shape,dtype=int)
        Z = zLs.sum(); Nseg = N-len(zLs)     # N nodes could only define N-1-(len(zLs)-1)=N-len(zLs) segments
        NsegRemain = Nseg - len(zLs)
        for i in range(len(Ns)):
            Ns[i] = zLs[i]/Z*NsegRemain
            Z -= zLs[i]; NsegRemain -= Ns[i]
        Ns = Ns+2
        return Ns
    while N < 2*cptmat.shape[1]:
        N *= 2
    Ns = calNs(zLs,N)

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

cvcpt = cpt2cmap(os.path.dirname(__file__)+'/cpts/cv_original.cpt','cv')[0]
rbcpt = cpt2cmap(os.path.dirname(__file__)+'/cpts/rainbow.cpt','rainbow')[0].reversed()

def buildBoundaryNorm(sep,N=256,n=None):
    n = N//(len(sep)-1) if n is None else n

    boundaries = np.zeros(1+n*(len(sep)-1))
    boundaries[0] = sep[0]
    for i in range(len(sep)-1):
        boundaries[(i*n)+1:(i+1)*n+1] = np.linspace(sep[i],sep[i+1],n+1)[1:]
    return mpl.colors.BoundaryNorm(boundaries,N)

def cmapZoom(cmap,zOld=[0,0.45,0.55,1],zNew=[0,0.25,0.75,1],N=1024):
    from matplotlib import cm
    from matplotlib.colors import ListedColormap
    samples = []
    for i in range(len(zOld)-1):
        samples += list(np.linspace(zOld[i],zOld[i+1],int(round((zNew[i+1]-zNew[i])*N))))
    cmap = cm.get_cmap(cmap,N*2) if type(cmap) is str else cmap
    return ListedColormap(cmap(samples))


''' miscellaneous '''

def colorCycle():
    fig = plt.figure()
    ax = fig.add_subplot()
    cycle = ax._get_lines.prop_cycler
    cl = [next(cycle) for _ in range(10)]
    cl = [c['color'] for c in cl]
    plt.close()
    return cl

def doubleArrow(x1,y1,x2,y2,**kwargs):
    plt.arrow(x1,y1,x2-x1,y2-y1,**kwargs)
    plt.arrow(x2,y2,x1-x2,y1-y2,**kwargs)


''' for plot objects in package shapely '''
def plot_coords(ax, ob, color=None):
    x, y = ob.xy
    ax.plot(x, y, 'o', color=color, zorder=1)

def plot_line(ax, ob, color=None, zorder=1, linewidth=3, alpha=1):
    x, y = ob.xy
    ax.plot(x, y, color=color, linewidth=linewidth, solid_capstyle='round', zorder=zorder, alpha=alpha)







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
