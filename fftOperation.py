import pyfftw
import numpy as np


'''         FFTW:                           Y2F or F2Y
if N=2l+1   f: 0 1 ... l -l ... -1          f: -l ... -1 0 1 ... l
if N=2l+2   f: 0 1 ... l l+1 -l ... -1      f: -l ... -1 0 1 ... l l+1
'''

def Y2F(dt,Y,norm=False,fftw=False): 
    ''' 
    In fact, normalization is meanless. I just make FFT of Asin(wt) to be A at w.
    '''
    N = len(Y)
    fNyq = 1/(2*dt)
    df = 1/(N*dt)
    Fraw = pyfftw.interfaces.numpy_fft.fft(Y) if fftw else np.fft.fft(Y)
    fraw = np.arange(N)*df
    F,f = np.zeros(N,dtype=complex),np.zeros(N)
    F[-(N//2+1):] = Fraw[:N//2+1]       # possitive  f branch
    F[:-(N//2+1)] = Fraw[N//2+1:]       # negative f branch
    f[-(N//2+1):] = fraw[:N//2+1]       
    f[:-(N//2+1)] = fraw[N//2+1:] - 2*fNyq
    if norm:
        F[f==0] = F[f==0]/N
        F[f!=0] = F[f!=0]/N*2
        F[abs(f)==fNyq] = F[abs(f)==fNyq]/2
    return f,F

def F2Y(f,F,fftw=False):
    N = len(F)
    try:
        df = f[1]-f[0]
    except:
        df = f
    dt = 1/(N*df)
    # fNyq = 1/(2*dt)
    Fraw,fraw = np.zeros(N,dtype=complex),np.zeros(N)
    Fraw[:N//2+1] = F[-(N//2+1):]
    Fraw[N//2+1:] = F[:-(N//2+1)]
    # fraw[:N//2+1] = f[-(N//2+1):]
    # fraw[N//2+1:] = f[:-(N//2+1)] + 2*fNyq
    Y = pyfftw.interfaces.numpy_fft.ifft(Fraw) if fftw else np.fft.ifft(Fraw)
    return dt,Y

def fillNegtiveWing(F):
    l = (len(F)-1)//2
    F.real[:l] = F.real[2*l:l:-1]
    F.imag[:l] = -F.imag[2*l:l:-1]
    return F

def inte(dt,Y,water_level=None):
    N = len(Y)
    L = N*dt
    f,F = Y2F(dt,Y)
    if water_level is not None:
        wl = water_level
        fwl = abs(f).max()/10**(wl/10)
        tmp = np.array([f0 if abs(f0)>fwl else fwl*np.sign(f0) for f0 in f[f!=0]])
        F[f!=0] = F[f!=0]/(2*1j*np.pi*tmp)
    else:
        F[f!=0] = F[f!=0]/(2*1j*np.pi*f[f!=0])
    F[f==0] = 0
    if N%2 == 0:
        F[-1] = 0
    _,DY = F2Y(f,F)
    return DY.real
    
def narrowBandGaussFilter(dt,Y,freq,alpha=20):
    f,F = Y2F(dt,Y)
    F = F * np.exp(-alpha*(f-freq)**2/freq**2)
    _,nBY = F2Y(f,F)
    return nBY

def randPhase(f,spectrum):
    F = spectrum.astype('complex')
    if len(F) % 2 == 0:
        raise Exception('We have not prepared for even F')
    arg = np.random.random(len(F[f>0]))*2*np.pi
    a = abs(F[f>0])*np.cos(arg)
    b = abs(F[f>0])*np.sin(arg)
    F[f>0] = a+b*1j
    F[f<0] = a[::-1]-b[::-1]*1j
    return F


