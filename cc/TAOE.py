#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 09:34:41 2020

@author: lun
"""

from obspy import UTCDateTime
from obspy.taup import TauPyModel
from obspy.clients.fdsn import Client
import obspy, os
import distaz, geo
import numpy as np

model = TauPyModel(model="iasp91")
fdsn_client = Client('IRIS')

evtfile = '../plot_doc/event_h6.lst'
evtorigin_lst = np.loadtxt(evtfile,usecols=2,dtype='str').tolist()
evt_lonlatdep =  np.loadtxt(evtfile,usecols=(3,4,5))
maglst = np.loadtxt(evtfile,usecols=1)
stlo = -140.148
stla = -8.855
#dur = 1800
time_beforetauP = 70
time_aftertauP = 30 

for i in range(len(evtorigin_lst)):
    evtorigin = evtorigin_lst[i]
    t1 = UTCDateTime(evtorigin)
    sacname = 'TAOE/%04d%02d%02d%02d%02d%02d_TAOE.BHZ'%(t1.year,t1.month,t1.day,t1.hour,t1.minute,t1.second)
    result = distaz.DistAz(stla,stlo,evt_lonlatdep[i,1],evt_lonlatdep[i,0])
    gcarc = result.getDelta()
    az = result.getAz()
    baz = result.getBaz()
    dist = result.degreesToKilometers(gcarc)
    try:
        arrivals = model.get_travel_times(source_depth_in_km=evt_lonlatdep[i,2],distance_in_degree=gcarc,phase_list=['P'])
        t_P = arrivals[0].time
        rayp = geo.srad2skm(arrivals[0].ray_param)
        st = fdsn_client.get_waveforms(network='G', station='TAOE', location='00',
                               channel='BHZ', starttime=t1+t_P-time_beforetauP, endtime=t1+t_P+time_aftertauP,
                               attach_response=True)
    except:
        print('%s: no P arrival'%sacname)
        continue
    pre_filt = (0.005, 0.01, 20.0, 25.0)
    st.remove_response(output='VEL', pre_filt=pre_filt)
    st.detrend(type='constant')
    st.detrend(type='linear')
    st.write(sacname,'SAC')
    try:
        st = obspy.read(sacname)[0]
        st.stats.sac.evdp = evt_lonlatdep[i,2]
        st.stats.sac.evla = evt_lonlatdep[i,1]
        st.stats.sac.evlo = evt_lonlatdep[i,0]
        st.stats.sac.stla = stla
        st.stats.sac.stlo = stlo
        st.stats.sac.mag = maglst[i]
        st.stats.sac.gcarc = gcarc
        st.stats.sac.dist = dist
        st.stats.sac.az = az
        st.stats.sac.baz = baz
        st.stats.sac.user0 = t1
        st.stats.sac.user1 = time_beforetauP
        st.stats.sac.user2 = time_aftertauP
        st.stats.sac.user3 = t_P
        st.stats.sac.user4 = rayp
        st.write(sacname,'SAC')
    except:
        print("can't read "+sacname)
        os.system('rm '+sacname)
        continue
    print(sacname)