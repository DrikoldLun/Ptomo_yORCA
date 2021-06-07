#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 19:41:15 2020

@author: lun
"""

import numpy as np
import os, obspy, configparser, sys
config = configparser.RawConfigParser()
import cc
    
def preprocess(tr,be=[-70,30],freq=[0.1,2],win=[-3,5]):
    npts = tr.stats.npts
    time = np.linspace(be[0],be[1],num=npts,endpoint=True)
    tr.filter('bandpass',freqmin=freq[0],freqmax=freq[1],zerophase=True)
    tr.data = tr.data[np.where((time>=win[0])&(time<=win[1]))]
    tr.taper(max_percentage=0.05,type='hann')
    return tr.data


def remccc(comp,mag):
    head = '../paraP.cfg'
    if os.path.isfile(head):
        config.read(head)
    else:
        print("head file doesn't exist!")
    prefix = head.replace(os.path.basename(head),'')
    out_path = os.path.join(prefix,config.get('path', 'out_path'))
    path = {}
    path['P'] = os.path.join(out_path,config.get('path', 'P_dir'))
    path['Z'] = os.path.join(out_path,config.get('path', 'Z_dir'))
    time_before = config.getfloat('para', 'time_beforetauP')
    time_after = config.getfloat('para', 'time_aftertauP')
    
    Zlst = np.loadtxt('../evtlog/evt_Z.log',usecols=0,dtype='str')
    Zmaglog = np.loadtxt('../evtlog/evt_Z.log',usecols=4)
    Zlst = Zlst[Zmaglog>=6]
    dat = 'cclog/%s_mag_%.1f_%.1f/%s_mag_%.1f_%.1f.dat'%(comp,mag[0],mag[1],comp,mag[0],mag[1])
    wf_file = 'cclog/%s_mag_%.1f_%.1f/%s_mag_%.1f_%.1f_winfreq.dat'%(comp,mag[0],mag[1],comp,mag[0],mag[1])
    wf_file_Z = 'cclog/Z_mag_6.0_10.0/Z_mag_6.0_10.0_winfreq.dat'
    #Pdat = '/home/lun/Desktop/Ptomo/cc/cclog/P_mag_6.0_10.0/P_mag_6.0_10.0.dat'
    #Pwf_file = '/home/lun/Desktop/Ptomo/cc/cclog/P_mag_6.0_10.0/P_mag_6.0_10.0_winfreq.dat'
    Cevtid = np.loadtxt(dat,usecols=0,dtype='str')
    wf = np.loadtxt(wf_file)
    wf_Z = np.loadtxt(wf_file_Z)

    lines = []
    with open(dat,'r') as f:
        Clines = f.readlines()
        f.close()
    Clines = np.array(Clines)

    evtid_uni = np.unique(Cevtid)
    #common_lst = [tmp for tmp in Zevtid_uni if tmp in Pevtid_uni]
    stalst = np.loadtxt(dat,usecols=4,dtype='str')
    lst = np.loadtxt('../evtlog/evt_%s.log'%comp,usecols=0,dtype='str')
    maglog = np.loadtxt('../evtlog/evt_%s.log'%comp,usecols=4)
    lst = lst[(maglog>=mag[0])&(maglog<mag[1])]
    #Zdtlst = np.loadtxt(dat,usecols=10)
    #Zccmaxlst = np.loadtxt(dat,usecols=-3)

    #lun ievtZ = np.where(Zlst==evtid)[0][0]
    
    for evtid in evtid_uni:
        #print(evtid)
        ievt = np.where(lst==evtid)[0][0]
        if evtid in Zlst:
            ievtZ = np.where(Zlst==evtid)[0][0]
            print(evtid)
        else:
            ievtZ = None
        #print(Zwf[ievt])
        Cstalst = stalst[Cevtid==evtid]
        if len(Cstalst) <= 4: continue
        Cline = Clines[Cevtid==evtid]
        data = []
        for sta in Cstalst:
            tr = obspy.read(os.path.join(path[comp],evtid,evtid+'_'+sta+'.%s'%comp))[0]
            if ievtZ == None:
                data.append(preprocess(tr,be=[-time_before,time_after],freq=wf[ievt,2:],win=wf[ievt,:2]))
            else:
                data.append(preprocess(tr,be=[-time_before,time_after],freq=wf[ievt,2:],win=wf_Z[ievtZ,:2]))
        data = np.array(data).T
        dcor,dcstd,dvcstd,ccmax,rsstd = cc.mccc(data,1./tr.stats.sampling_rate,ifwt=True)
        #data = matlab.double(data.tolist())
        #[dcor,dcstd,dvcstd,resid] = eng.xcortimes_ze(data,1./tr.stats.sampling_rate,'lagmax',1.,'ifwt',bool(0),nargout=4)
        for i in range(len(Cstalst)):
            line = Cline[i].split()
            newline = '  '.join(line[:10]+['%.3f'%dcor[i],'%.4f'%dcstd[i],'%.4f'%dvcstd[i],'%.4f'%ccmax[i],'%.5f'%rsstd[i]]+line[15:])+'\n'
            lines.append(newline)
        
        with open('%s_mag_%.1f_%.1f.dat'%(comp,mag[0],mag[1]),'w+') as f:
            f.writelines(lines)
            f.close()

def main():
    head = '../paraP.cfg'
    if os.path.isfile(head):
        config.read(head)
    else:
        print("head file doesn't exist!")
    prefix = head.replace(os.path.basename(head),'')
    out_path = os.path.join(prefix,config.get('path', 'out_path'))
    path = {}
    path['P'] = os.path.join(out_path,config.get('path', 'P_dir'))
    path['Z'] = os.path.join(out_path,config.get('path', 'Z_dir'))
    time_before = config.getfloat('para', 'time_beforetauP')
    time_after = config.getfloat('para', 'time_aftertauP')
    Zdat = 'cclog/Z_mag_6.0_10.0/Z_mag_6.0_10.0.dat'
    Zwf_file = 'cclog/Z_mag_6.0_10.0/Z_mag_6.0_10.0_winfreq.dat'
    Pdat = 'cclog/P_mag_6.0_10.0/P_mag_6.0_10.0.dat'
    Pwf_file = 'cclog/P_mag_6.0_10.0/P_mag_6.0_10.0_winfreq.dat'
    Zevtid = np.loadtxt(Zdat,usecols=0,dtype='str')
    Zwf = np.loadtxt(Zwf_file)
    Pevtid = np.loadtxt(Pdat,usecols=0,dtype='str')
    Pwf = np.loadtxt(Pwf_file)

    lines = []
    with open(Zdat,'r') as f:
        Zlines = np.array(f.readlines())
        f.close()
    with open(Pdat,'r') as f:
        Plines = np.array(f.readlines())
        f.close()

    Zevtid_uni = np.unique(Zevtid)
    Pevtid_uni = np.unique(Pevtid)
    common_lst = [tmp for tmp in Zevtid_uni if tmp in Pevtid_uni]
    Zstalst = np.loadtxt(Zdat,usecols=4,dtype='str')
    Pstalst = np.loadtxt(Pdat,usecols=4,dtype='str')
    Zlst = np.loadtxt('../evtlog/evt_Z.log',usecols=0,dtype='str')
    Zmaglog = np.loadtxt('../evtlog/evt_Z.log',usecols=4)
    Zlst = Zlst[Zmaglog>=6]
    #Zdtlst = np.loadtxt(Zdat,usecols=10)
    #Zccmaxlst = np.loadtxt(Zdat,usecols=-3)
    Plst = np.loadtxt('../evtlog/evt_P.log',usecols=0,dtype='str')
    Pmaglog = np.loadtxt('../evtlog/evt_P.log',usecols=4)
    Plst = Plst[Pmaglog>=6]
    #Pdtlst = np.loadtxt(Pdat,usecols=10)
    
    # common events
    for evtid in common_lst:
        #if evtid != '20180916211148':
        #    continue
        print(evtid)
        ievtZ = np.where(Zlst==evtid)[0][0]
        staZ = Zstalst[Zevtid==evtid].tolist()
        Zdata = []
        for sta in staZ:
            tr = obspy.read(os.path.join(path['Z'],evtid,evtid+'_'+sta+'.Z'))[0]
            Zdata.append(preprocess(tr,be=[-time_before,time_after],freq=Zwf[ievtZ,2:],win=Zwf[ievtZ,:2]))
        Zdata = np.array(Zdata).T
        dtZ = 1./tr.stats.sampling_rate
        ievtP = np.where(Plst==evtid)[0][0]
        staP = Pstalst[Pevtid==evtid].tolist()
        Pdata = []
        for sta in staP:
            tr = obspy.read(os.path.join(path['P'],evtid,evtid+'_'+sta+'.P'))[0]
            Pdata.append(preprocess(tr,be=[-time_before,time_after],freq=Pwf[ievtP,2:],win=Zwf[ievtZ,:2]))
        Pdata = np.array(Pdata).T
        dtP = 1./tr.stats.sampling_rate
        try:
            sta,dcor,dcstd,dvcstd,rsstd = cc.ZP_mccc(staZ,Zdata,dtZ,staP,Pdata,dtP,lagmax=1.,ifwt=True)
        except:
            print('singular matrix: continue')
            continue
        Zline = Zlines[Zevtid==evtid]
        Pline = Plines[Pevtid==evtid]
        sta,staZ,staP = map(np.array,[sta,staZ,staP])
        for i in range(len(sta)):
            if sta[i] in staZ:
                line = Zline[staZ==sta[i]][0]
            else:
                line = Pline[staP==sta[i]][0]
            line = line.split()
            newline = '  '.join(line[:10]+['%.3f'%dcor[i],'%.4f'%dcstd[i],'%.4f'%dvcstd[i],'0','%.5f'%rsstd[i],line[-2],'Z&P'])+'\n'
            lines.append(newline)
    # vertical events
    for evtid in Zevtid_uni:
        if evtid in common_lst: continue
        lines += Zlines[Zevtid==evtid].tolist()
    # pressure events       
    for evtid in Pevtid_uni:
        if evtid in common_lst: continue
        lines += Plines[Pevtid==evtid].tolist()
    # extra traces from low magnitude    
    Zdat_add = 'Z_mag_5.5_6.0.dat'
    with open(Zdat_add,'r') as f:
        lines += f.readlines()
        f.close()
    
    PZlst = np.loadtxt('../evtlog/evt_PorZ.log',usecols=0,dtype='str')
    newlines = []
    for line in lines:
        line = line.split()
        newline = '  '.join(line[:6]+['%d'%np.where(PZlst==line[0])[0]]+line[7:])+'\n'
        newlines.append(newline)
    with open('PZmix.dat','w+') as f:
        f.writelines(newlines)
        f.close()
    
if __name__ == '__main__':
    #remccc('Z',[5.5,6])
    main()
