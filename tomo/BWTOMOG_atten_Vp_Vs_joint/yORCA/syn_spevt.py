#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 08:41:26 2020

@author: lun
"""
import random
import numpy as np
import distaz
from obspy.taup import TauPyModel
import geo

def fake_event(total,center,distrange):
    # center [lat lon]
    # dist [deg]
    evt_info = []
    for ie = range(length(total)):
        azi = ie%360
        dist = random.uniform(distrange[0],distrange[1])
        [lat, lon] = latlon_from(center[0], center[1], azi, dist)
        evt_info.append()


def sample_seg(total,minimum,maximum):
    # total, minimum, maximum should be interger
    s = []
    for i in range(total):
        if len(s) > 0 and s[-1] > maximum:
            s = s[:-1]
        elif sum(s) == total:
            break
        elif sum(s)+minimum >= total:
            s = s[:-1]
            s.append(total-sum(s))
        else:
            s.append(random.randint(minimum,maximum))
    return s

total = 8000
cf = 0.3
evtfile = '../../../evtlog/evt_Z.log'
stafile = 'data/stations_PZ.dat'

stantw = np.loadtxt(stafile,usecols=(0,1),delimiter=',',dtype='str') #sta, ntw
stainfo = np.loadtxt(stafile,usecols=(2,3,4),delimiter=',') # stla, stlo, elev
evtlog = np.loadtxt(evtfile,usecols=(1,2,3)) #evlo, evla, evdp
stanum = stantw.shape[0]
evtnum = evtlog.shape[0]

trnum_lst = sample_seg(total,4,12) # traces for each event
evtid = random.sample(range(evtnum),len(trnum_lst)) # sample events
evtlog = evtlog[evtid]
lines = []
model = TauPyModel(model="iasp91")

for i in range(len(evtid)):
    evlo, evla, evdp = evtlog[i]
    staid = random.sample(range(stanum),trnum_lst[i]) # sample stations for each event
    for j in staid:
        stla, stlo, elev = stainfo[j]
        result = distaz.DistAz(stla,stlo,evla,evlo)
        baz = result.getBaz()
        gcarc = result.getDelta()
        try:
            arrivals = model.get_travel_times(source_depth_in_km=evdp,distance_in_degree=gcarc,phase_list=['P'])
            rayp = np.deg2rad(geo.srad2skm(arrivals[0].ray_param)*6371)
            line = '%8.4f %8.3f %8.2f %8s %8s %8d %8.4f %8.4f %8.2f %8.3f %8.4f %8.4f %8.4f\n'%(rayp,gcarc,baz,stantw[j,0],stantw[j,1],i,evla,evlo,evdp,0,0,1,cf)
            lines.append(line)
        except:
            continue

synevtfile = 'data/syn%d_dT_PZ.dat'%len(lines)
print(len(lines))

with open(synevtfile,'w+') as f:
    f.writelines(lines)
    f.close()
