#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 01:47:49 2020

@author: lun
"""
import numpy as np
#datdir = 'cclog'
#suffix = '_mag_6.0_10.0'
#Pdat = os.path.join(datdir,'P'+suffix,'P'+suffix+'.dat')
#Zdat = os.path.join(datdir,'Z'+suffix,'Z'+suffix+'.dat')
mode = 'direct' #direct(only P or Z measurements for 1 evt),mix
dominant = 'both' # Z,P,num,both

Zdat = 'Z_mag_6.0_10.0.dat'
Pdat = 'P_mag_6.0_10.0.dat'

Zevtid = np.loadtxt(Zdat,usecols=0,dtype='str')
Pevtid = np.loadtxt(Pdat,usecols=0,dtype='str')
Zevtid_uni = np.unique(Zevtid)
Pevtid_uni = np.unique(Pevtid)
PZlst = np.loadtxt('../evtlog/evt_PorZ.log',usecols=0,dtype='str')
common_lst = [tmp for tmp in Zevtid_uni if tmp in Pevtid_uni]
Zstalst = np.loadtxt(Zdat,usecols=4,dtype='str')
Pstalst = np.loadtxt(Pdat,usecols=4,dtype='str')
Zdtlst = np.loadtxt(Zdat,usecols=10)
Pdtlst = np.loadtxt(Pdat,usecols=10)

lines = []

with open(Zdat,'r') as f:
    Zlines = f.readlines()
    f.close()
with open(Pdat,'r') as f:
    Plines = f.readlines()
    f.close()

for line in Plines:
    if line.split()[0] not in common_lst:
        lines.append(line)
for line in Zlines:
    if line.split()[0] not in common_lst:
        lines.append(line)

if mode == 'direct':
    for evtid in common_lst:
        Zid = np.where(Zevtid==evtid)[0]
        Pid = np.where(Pevtid==evtid)[0]
        if dominant == 'num':
            if len(Zid) >= len(Pid):
                lines += [Zlines[i] for i in Zid]
            else:
                lines += [Plines[i] for i in Pid]
        elif dominant == 'Z':
            lines += [Zlines[i] for i in Zid]
        elif dominant == 'P':
            lines += [Plines[i] for i in Pid]
        elif dominant == 'both':
            lines += [Zlines[i] for i in Zid]
            lines += [Plines[i] for i in Pid]

'''        
elif mode == 'mix':
    for evtid in common_lst:
        Zid = np.where(Zevtid==evtid)[0]
        Zsta = Zstalst[Zid].tolist()
        #Zdt = Zdtlst[Zid]
        Pid = np.where(Pevtid==evtid)[0]
        Psta = Pstalst[Pid].tolist()
        Pdt = Pdtlst[Pid]
        common_sta = [sta for sta in Zsta if sta in Psta]
        if common_sta == []:
            if dominant == 'num':
                if len(Zid) >= len(Pid):
                    lines += [Zlines[i] for i in Zid]
                else:
                    lines += [Plines[i] for i in Pid]
            elif dominant == 'Z':
                lines += [Zlines[i] for i in Zid]
            elif dominant == 'P':
                lines += [Plines[i] for i in Pid]
        else:
            Zid_common = np.array([Zid[Zsta.index(sta)] for sta in common_sta])
            Zid_rest = np.array([i for i in Zid if i not in Zid_common]).astype('int')
            Pid_common = np.array([Pid[Psta.index(sta)] for sta in common_sta])
            Pid_rest = np.array([i for i in Pid if i not in Pid_common]).astype('int')
            Zstatic = np.sum(Zdtlst[Zid_common])
            Pstatic = np.sum(Pdtlst[Pid_common])
            shift1 = (Zstatic-Pstatic)/len(common_sta)
            Pdt += shift1
            if dominant == 'Z':
                Zid_merge, Pid_merge = Zid[:], Pid_rest[:] 
            elif dominant == 'P':
                Zid_merge, Pid_merge = Zid_rest[:], Pid[:] 
            elif dominant == 'num':
                if len(Zid) >= len(Pid):
                    Zid_merge, Pid_merge = Zid[:], Pid_rest[:]
                else:
                    Zid_merge, Pid_merge = Zid_rest[:], Pid[:]
            shift2 = -np.mean(np.hstack((Zdtlst[Zid_merge],Pdtlst[Pid_merge])))
            print(Zdtlst[Zid])
            Zdtlst[Zid] = Zdtlst[Zid]+shift2
            Pdtlst[Pid] = Pdtlst[Pid]+shift2
            print(sum(Zdtlst[Zid_merge])+sum(Pdtlst[Pid_merge]))
            for i in Zid_merge:
                line = Zlines[i].split()
                lines.append('  '.join(line[:10]+['%.3f'%Zdtlst[i]]+line[11:]))
            for i in Pid_merge:
                line = Plines[i].split()
                lines.append('  '.join(line[:10]+['%.3f'%Pdtlst[i]]+line[11:]))
'''

Zdat_add = 'Z_mag_5.5_6.0.dat'
with open(Zdat_add,'r') as f:
    lines += f.readlines()
    f.close()

newlines = []
for line in lines:
    line = line.replace('\n','').split()
    if line[-1] == 'Z':
        newline = '  '.join(line[:6]+['%d'%np.where(PZlst==line[0])[0]]+line[7:])+'\n'
    elif line[-1] == 'P':
        newline = '  '.join(line[:6]+['-%d'%np.where(PZlst==line[0])[0]]+line[7:])+'\n'
    newlines.append(newline)
    
with open('PZ_direct_%s.dat'%dominant,'w+') as f:
   f.writelines(newlines)
   f.close()
