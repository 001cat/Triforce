import obspy,time
from Triforce.pltHead import *
from Triforce.fftOperation import Y2F
from Triforce.obspyPlus import rmRESPYe

''' some pitfalls '''
# st = obspy.read('Data/CCM.IU.mseed')
# st.attach_response(obspy.read_inventory('Data/IRISDMC-CCM.IU.xml'))

# st.trim(obspy.UTCDateTime('20110801'),obspy.UTCDateTime('20110803'))
# st.detrend('linear')
# st.merge(fill_value=0)

# tr0 = st[0].copy()  # Ye's method
# tr1 = st[0].copy()  # similar result from obspy, need pre_filt and no taper
# tr2 = st[0].copy()  # wrong result from without pre_filt

# tr0 = rmRESPYe(tr0,'Data/IRISDMC-CCM.IU.xml')
# f2,f3 = 0.01*0.7,0.2*1.3
# f1,f4 = f2*0.8,f3*1.2
# tr1.remove_response(pre_filt=[f1,f2,f3,f4],taper=False)
# tr2.remove_response(taper=False)
# plt.figure()
# for tr in [tr0,tr1,tr2]:
#     plt.plot(tr.times(),tr.data)

# # To apply filter, move the first point back to zero is necessary
# tr3 = tr2.copy()
# tr2.data -= tr2.data[0];tr2.filter('bandpass',freqmin=0.01,freqmax=0.2,zerophase=True) 
# tr3.filter('bandpass',freqmin=0.01,freqmax=0.2,zerophase=True)
# plt.figure()
# for tr in [tr0,tr2,tr3]:
#     plt.plot(tr.times(),tr.data)
# # not only the initial tremble, but also the remaining parts are different
# # the spectrum will also become totally different


''' Speed test '''
st = obspy.read('Data/CCM.IU.mseed')
st.attach_response(obspy.read_inventory('Data/IRISDMC-CCM.IU.xml'))

# st.trim(obspy.UTCDateTime('20110801'),obspy.UTCDateTime('20110803')) # short time events
st.trim(obspy.UTCDateTime('20110720'),obspy.UTCDateTime('20111008')) # long time series

st.detrend('linear')
st.merge(fill_value=0)

tr0 = st[0].copy()
tr1 = st[0].copy() 
tr2 = st[0].copy()
tr3 = st[0].copy()

f2,f3 = 0.01*0.7,0.2*1.3
f1,f4 = f2*0.8,f3*1.2
t0 = time.time(); tr0.remove_response(pre_filt=[f1,f2,f3,f4],taper=False); print(time.time()-t0)
t0 = time.time(); tr1.remove_response(taper=False); tr1.data -= tr1.data[0]; print(time.time()-t0)
t0 = time.time(); tr2 = rmRESPYe(tr2,'Data/IRISDMC-CCM.IU.xml'); print(time.time()-t0)
t0 = time.time(); tr3 = rmRESPYe(tr3,'Data/IRISDMC-CCM.IU.xml',fftw=True); print(time.time()-t0)

# tr1.filter('bandpass',freqmin=0.01,freqmax=0.2,zerophase=True) 
# plt.figure()
# for tr in [tr0,tr1,tr2,tr3]:
#     plt.plot(tr.times(),tr.data)

stDeNoised = obspy.Stream([tr0,tr1,tr2,tr3])
stDeNoised.filter('bandpass',freqmin=0.01,freqmax=0.2,zerophase=True) 
plt.figure()
for tr in stDeNoised:
    plt.plot(tr.times(),tr.data)

''' 
short 2*86400 pts 
0.4328036308288574
0.38201117515563965
0.5798664093017578
0.6069834232330322

long 80*86400 pts
27.69945216178894
27.702144861221313
11.46040964126587
6.647719383239746
'''
