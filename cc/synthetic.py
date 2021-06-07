#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  4 23:48:56 2020

@author: lun
"""
import obspy, os
from obspy import UTCDateTime
from obspy.taup import TauPyModel
import numpy as np
import distaz, geo

model = TauPyModel(model="iasp91")
evtfile = '../plot_doc/event_h6.lst'
maglst = np.loadtxt(evtfile,usecols=1)
evtorigin_lst = np.loadtxt(evtfile,usecols=2,dtype='str').tolist()
evt_lonlatdep =  np.loadtxt(evtfile,usecols=(3,4,5))
time_beforetauP = 70
time_aftertauP = 30

stla, stlo = np.loadtxt('rec1.txt',usecols=(0,1)) 

if not os.path.exists('synthetic'):
    os.mkdir('synthetic')

for i in range(len(evt_lonlatdep)):
    evtorigin = evtorigin_lst[i]
    t1 = UTCDateTime(evtorigin)
    ymdhs = '%4d%02d%02d%02d%02d'%(t1.year,t1.month,t1.day,t1.hour,t1.minute)
    #ymdhs = '201809061549'
    evttime = '%s%02d'%(ymdhs,t1.second)
    sacname = 'synthetic/%s_syn.BHZ'%evttime
    
    result = distaz.DistAz(stla,stlo,evt_lonlatdep[i,1],evt_lonlatdep[i,0])
    gcarc = result.getDelta()
    az = result.getAz()
    baz = result.getBaz()
    dist = result.degreesToKilometers(gcarc)
    try:
        arrivals = model.get_travel_times(source_depth_in_km=evt_lonlatdep[i,2],distance_in_degree=gcarc,phase_list=['P'])
        t_P = arrivals[0].time
        rayp = geo.srad2skm(arrivals[0].ray_param)
    except:
        print('%s: no P arrival'%sacname)
        continue
    try:
        os.system('./FetchSyn -recfile rec1.txt -model iasp91_2s -units velocity -evid GCMT:%sA -sdepm %f'%(ymdhs,evt_lonlatdep[i,2]*1000))
        os.system('unzip Synthetics.zip')
        os.system('mv *Z.sac %s'%sacname)
        os.system('rm Synthetics.zip Syngine.log *.sac')
    except:
        print("%s: can't proceed downloading")
        continue
    try:
        st = obspy.read(sacname)[0]
        nb, ne = int((t_P-time_beforetauP)*st.stats.sampling_rate), int((t_P+time_aftertauP)*st.stats.sampling_rate)+1
        st.data = st.data[nb:ne]
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