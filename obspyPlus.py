import obspy,os
import numpy as np
import matplotlib.pyplot as plt
from Triforce.fftOperation import Y2F,F2Y
from Triforce.utils import get_current_memory,randString

def midnightBefore(t):
    return obspy.UTCDateTime(t.year,t.month,t.day)

def genDateTimeWins(starttime,endtime,winL=86400):
    starttime = obspy.UTCDateTime(starttime)
    endtime   = obspy.UTCDateTime(endtime)
    t0 = starttime
    winLst = []
    while t0< endtime:
        winLst.append((t0,t0+winL))
        t0 += winL
    return winLst

def seismicPreProc(stIn,inv,starttime=None,endtime=None,freqmin=0.01,freqmax=None, 
                   reSampleDelta=False,merge=True,rmResp=True,zerophase=True,**kwargs):
    ''' Stream pre-processing including: cut, merge, resample, filter, remove response

    Args:
        stIn,inv:           input stream and inventory, original stream won't be changed
        starttime,endtime:  start and end time of output
        reSampleDelta:      delta time of output, an anti-alias lowpass filter would be applied
        rmResp:             remove response
        freqmin,freqmax:    filter parameters
        merge:              merge traces with same seed ID, set as 'highpass' to apply 0.01Hz 
                            highpass filter after merge
    Other procedure:
        detrend:            default is 'linear', set detrend='simple' to avoid out of memory
        taper:              valid when filter would be applied
    '''
    st = stIn.copy()

    if reSampleDelta or rmResp or (freqmin is not None) or (freqmax is not None):
        border = 500 if 'border' not in kwargs.keys() else kwargs['border']
    else:
        border = 0
    
    for i in range(len(st)-1,-1,-1):
        if st[i].stats.npts < 2*border:
            del st[i]

    starttime = stStartEndTime(st,True)[0]+border if starttime is None else obspy.UTCDateTime(starttime)
    endtime   = stStartEndTime(st,True)[1]-border if endtime is None else obspy.UTCDateTime(endtime)
    st.trim(starttime-border,endtime+border)

    try:
        st.detrend(type=kwargs['detrend'])
    except:
        st.detrend(type='linear')
    st.taper(None,max_length=border)
    if reSampleDelta:
        # print('reSampling ...')
        st.filter('lowpass',corners=6,freq=0.5/reSampleDelta/2,zerophase=True)
        for tr in st:
            if (starttime-border) >= tr.stats.starttime:
                interpStart = starttime-border
            else:
                interpStart = starttime-border + np.ceil((tr.stats.starttime-(starttime-border))/reSampleDelta)*reSampleDelta
            tr.interpolate(1.0/reSampleDelta,starttime=interpStart)
    if rmResp:
        # print('removing response ...')
        st.attach_response(inv)
        for tr in st:
            if tr.stats.response.response_stages[0].input_units == 'PA':
                outputUnit = 'VEL'
            elif tr.stats.response.response_stages[0].input_units in ('M/S','M'):
                outputUnit='DISP'
            else:
                print(f'Unknown input unit, no response removed: {tr.get_id()}')
                continue
            try:
                respYe = kwargs['respYe']
            except:
                respYe = False
            if respYe:
                tr = rmRESPYe(tr,inv,freqmin=0.002,freqmax=min(1000,1/tr.stats.delta/2),fftw=True)
            else:
                tr.remove_response(pre_filt=[0.001,0.002,1000,2000],output=outputUnit)

    if (freqmin is None) and freqmax:
        st.filter('lowpass',freq=freqmax,zerophase=zerophase)
    elif freqmin and (freqmax is None):
        st.filter('highpass',freq=freqmin,zerophase=zerophase)
    elif freqmin and freqmax:
        st.filter('bandpass',freqmin=freqmin,freqmax=freqmax,zerophase=zerophase)
    else:
        pass
 
    for i,tr in zip(range(len(st)-1,-1,-1),st[::-1]):
        minL = 2*border
        # if freqmin is not None:
        #     minL = max(minL,10/freqmin)
        if tr.stats.endtime - tr.stats.starttime > minL:
            tr.trim(tr.stats.starttime+border,tr.stats.endtime-border)
            pass
        else:
            del st[i]
            
    if merge is True:
        st.merge(fill_value='interpolate')
    elif merge == 'highpass':
        st.merge(fill_value='interpolate')
        st.filter('highpass',freq=0.01,zerophase=zerophase)
    return st

def rmRESPYe(trIn,xmlResp,freqmin=0.01,freqmax=0.2,copy=False,fftw=False):
    ''' remove response following Ye's code '''
    def flatHann(f,f1,f2,f3,f4):
        win = np.zeros(f.shape)
        I = (f>f1)*(f<=f2)
        win[I] = ((1-np.cos(np.pi*(f1-f)/(f2-f1))) * 0.5)[I]
        I = (f>f2)*(f<=f3)
        win[I] = 1
        I = (f>f3)*(f<=f4)
        win[I] = ((1+np.cos(np.pi*(f3-f)/(f4-f3))) * 0.5)[I]
        return win
    def fillNega(F):
        N = (len(F)-1)//2
        F.real[:N] = F.real[2*N:N:-1]
        F.imag[:N] = -F.imag[2*N:N:-1]
        return F
    tr = trIn.copy() if copy else trIn
    net,sta,loc,cha = tr.id.split('.')
    year,jday = tr.stats.starttime.year,tr.stats.starttime.julday
    f2,f3 = freqmin*0.7,freqmax*1.3
    f1,f4 = f2*0.8,f3*1.2
    xml2respBin = '/home/ayu/Packages/IRIS/xml2resp'
    evalrespBin = '/home/ayu/Packages/IRIS/evalresp'

    tmpDir = f'tmp{randString(10)}'
    os.system(f'mkdir {tmpDir}')
    if type(xmlResp) == obspy.Inventory:
        xmlResp.write(f'{tmpDir}/inv.xml',format='stationxml')
        xmlResp = f'{tmpDir}/inv.xml'
    os.system(f'{xml2respBin} -o {tmpDir}/RESP.{net}.{sta} {xmlResp} > /dev/null 2>&1')
    os.system(f'cd {tmpDir} && '+
    f'{evalrespBin} {sta} {cha} {year} {jday} {f1} {f4} 100 -f RESP.{net}.{sta} -v > /dev/null 2>&1')
    freq,amp = np.loadtxt(f'{tmpDir}/AMP.{tr.id}').T
    _,pha = np.loadtxt(f'{tmpDir}/PHASE.{tr.id}').T
    os.system(f'rm -r {tmpDir}')
    # amp *= 0.000000001
    pha *= np.pi/180

    tr.detrend(type='linear')
    f,F = Y2F(tr.stats.delta,tr.data,fftw=fftw)
    newF = np.zeros(F.shape,dtype=complex)
    I = (f>=f1) * (f<=f4)
    F = F[I]
    amp = np.interp(f[I], freq, amp)
    pha = np.interp(f[I], freq, pha)

    cosTaper = flatHann(f,f1,f2,f3,f4)

    newF.real[I] = (F.real*np.cos(pha)+F.imag*np.sin(pha))/amp
    newF.imag[I] = (F.imag*np.cos(pha)-F.real*np.sin(pha))/amp
    newF = newF * cosTaper
    F = fillNega(newF)
    _,Y = F2Y(f,F,fftw=fftw)
    tr.data = Y.real
    # tr.data = Y*2. / len(F)
    return tr



def stStartEndTime(stIn,secRound=False):
    '''Get first start time and last end time of stream'''
    def round2sec(t):
        return t-t.microsecond/1e6+round(t.microsecond/1e6)
    starttime = stIn[0].stats.starttime
    endtime   = stIn[0].stats.endtime
    for tr in stIn:
        if tr.stats.starttime < starttime:
            starttime = tr.stats.starttime
        if tr.stats.endtime > endtime:
            endtime = tr.stats.endtime
    if secRound:
        starttime,endtime = round2sec(starttime),round2sec(endtime)
    return starttime,endtime
def stStartEndTimeAllIn(stIn):
    '''Get last start time and first end time of stream'''
    starttime = stIn[0].stats.starttime
    endtime   = stIn[0].stats.endtime
    for tr in stIn:
        if tr.stats.starttime > starttime:
            starttime = tr.stats.starttime
        if tr.stats.endtime < endtime:
            endtime = tr.stats.endtime
    return starttime,endtime


def sepLargeMseed(mseedIn,seedID,outdir,starttime,endtime,segL=86400):
    ''' Seperate Large miniseed file into smaller ones that could be read by obspy 
    Args:
        mseedIn:            miniseed file path
        seedID:             the seed ID of extracted trace
        outdir:             output directory
        starttime,endtime:  only segments btween start and end time would be extacted
        segL:               segment length in seconds
        '''
    import io
    reclen = 512; chunksize = 250 * 2048 * reclen # 250 MB, depends on memory capacity 

    sID = seedID
    btime = obspy.UTCDateTime(starttime)
    etime = obspy.UTCDateTime(endtime)
    stetimePrevious = obspy.UTCDateTime(1900,1,1)
    with io.open(mseedIn, "rb") as fh:
        flag = False
        while True:
            with io.BytesIO() as buf:
                c = fh.read(chunksize)
                if not c:
                    break
                buf.write(c)
                buf.seek(0, 0)
                stRaw = obspy.read(buf)
                print(stRaw)
                st = stRaw.select(id=sID)
                if len(st) == 0:
                    continue
                stbtime, stetime = stStartEndTime(st)
                if stetime < btime:
                    continue
                while stetime >= btime + segL:
                    stN = st.copy()
                    stN.trim(btime,min(btime + segL,etime))
                    if len(stN) > 0:
                        stN.write(f'{outdir}/{sID}.{btime.strftime("%Y%m%d")}.mseed')
                        print(btime.strftime("%Y%m%d"))
                    btime = btime + segL
                    if btime >= etime:
                        break
                if stRaw[-1].id != sID or stetimePrevious==stetime:
                    stN = st.copy()
                    stN.trim(btime,min(btime + segL,etime))
                    stN.write(f'{outdir}/{sID}.{btime.strftime("%Y%m%d")}.mseed')
                    print(btime.strftime("%Y%m%d"))
                    btime = btime + segL
                    if stetimePrevious==stetime:
                        break
                if btime >= etime:
                        break
                stetimePrevious = stetime
                fh.seek(-chunksize//2,1)

def splitMseed(mseedIn,outDir=None,segMB=250):
    if outDir is None:
        outDir = '/'.join(mseedIn.split('/')[:-1])
    fname = '.'.join(mseedIn.split('/')[-1].split('.')[:-1])
    i = 0
    import io
    reclen = 512; chunksize = segMB * 2048 * reclen
    with io.open(mseedIn, "rb") as fh:
        while True:
            i = i+1
            c = fh.read(chunksize)
            if not c:
                break
            with io.open(f'{outDir}/{fname}-{i:03d}.mseed','wb') as buf:
                buf.write(c)
 


def sendDataRequest(ntstaList,starttime,endtime,comps,label='Test',
                    requestName='Ayu',email='ayu4cucssa@gmail.com'):
    ''' Request data through iris's BREQ_FAST service '''
    LABEL = label
    starttime = obspy.UTCDateTime(starttime)
    endtime = obspy.UTCDateTime(endtime)

    requestLines = []
    for ntsta in ntstaList:
        requestLines.append(f'{ntsta.split(".")[1]} {ntsta.split(".")[0]} '+
                            f'{starttime.strftime("%Y %m %d %H %M %S")} ' + 
                            f'{endtime.strftime("%Y %m %d %H %M %S")} {len(comps)} {" ".join(comps)}')

    with open(f'DataRequest.txt','w') as f:
        f.write(f'.Name {requestName}\n')
        f.write(f'.INST CU\n')
        f.write(f'.MAIL University of Colorado at Boulder\n')
        f.write(f'.EMAIL {email}\n')
        f.write(f'.PHONE\n.FAX\n.MEDIA: Electronic (FTP)\n.ALTERNATE MEDIA: Electronic (FTP)\n')
        f.write(f'.LABEL {LABEL}\n')
        f.write(f'.END\n')
        for line in requestLines:
            f.write(f'{line}\n')
    import os
    print('Sending mail')
    # if your isp allow port 25
    os.system(f"cat DataRequest.txt | mail -s 'Requesting Data' miniseed@iris.washington.edu")
    # if your isp do block port 25
    # os.system(f"echo 'To:miniseed@iris.washington.edu\nSubject: Data Request\n' "
    #           f" | cat - DataRequest.txt | msmtp -a ayu4cucssa -t")

def sendDataRequestByEvent(evt,inv,btime=0,etime=7200,comps=['LHZ'],
                           requestName='Ayu',email='ayu4cucssa@gmail.com'):
    t0 = evt.preferred_origin().time
    evtDate = evt.preferred_origin().time.strftime("%Y-%m-%d")
    evtMag  = f'{evt.preferred_magnitude().magnitude_type}{evt.preferred_magnitude().mag}'.lower().replace('.','')
    evtDescription = '-'.join(evt.event_descriptions[0].text.split()).lower()
    evtID = f'{evtDate}-{evtMag}-{evtDescription}'
    inv = inv.select(time=t0)
    ntstaList = [f'{net.code}.{sta.code}' for net in inv for sta in net]
    print(f'Sending data request: {evtID}')
    sendDataRequest(ntstaList,starttime=t0+btime,endtime=t0+etime,comps=comps,label=evtID,
                    requestName=requestName,email=email)



def STALTA(data,winL,winS):
    if (type(data) is not np.ndarray) or (np.any(np.isnan(data))):
        raise ValueError('Input data are invalid.')
    data = np.concatenate((np.ones(winL)*(np.sqrt(((data)**2).mean())),data))**2
    STA,LTA = np.zeros(data.shape),np.zeros(data.shape)
    for i in range(len(STA)):
        if i == 0:
            STA[i] = data[winL-winS+1:winL+1].mean()
            LTA[i] = data[1:winL+1].mean()
        else:
            STA[i] = STA[i-1] - data[winL-winS+i]/winS + data[winL+i]/winS
            LTA[i] = LTA[i-1] - data[i]/winL + data[winL+i]/winL
    return np.sqrt(STA/LTA)

def obspyFilter(dt,dataIn,type,**kwargs):
    tr = obspy.Trace()
    tr.stats.delta = dt
    tr.data = dataIn
    tr.detrend('linear')
    tr.taper(0.1)
    tr.filter(type,**kwargs)
    return tr.data




if __name__ == '__main__':
    plt.ion()
    st = obspy.read('/Users/ayu/Study/Data/LP26/try1-CCM.462971/CCM.IU.mseed')
    inv = obspy.read_inventory('/Users/ayu/Study/Data/LP26/try1-CCM.462971/IRISDMC-CCM.IU.xml')
    st = seismicPreProc(st,inv,'20040101T000000','20050101T000000')
    st.plot(handle=True)


