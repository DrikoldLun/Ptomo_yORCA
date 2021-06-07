import os, glob, obspy, configparser
import numpy as np


evtlog_path = 'evtlog'
if not os.path.exists(evtlog_path): os.mkdir(evtlog_path)

def evt_catalog(ch,ch_path,evtlog_path):
    # ch - P, Z, P&Z, P|Z
    evt_info = []
    if ch == 'P' or ch == 'Z':
        # evttime, evlo, evla, evdp, mag, trnum
        evt_lst = os.listdir(ch_path[ch])
        for evt in evt_lst:
            path_evt = os.path.join(ch_path[ch],evt)
            sac_lst = glob.glob(os.path.join(path_evt,evt+'*.'+ch))
            if len(sac_lst) == 0: os.rmdir(path_evt)
            else:
                stats = obspy.read(sac_lst[0])[0].stats.sac
                evlo, evla, evdp, mag = stats.evlo, stats.evla, stats.evdp, stats.mag
                stanum = len(sac_lst)
                evt_info.append([evt,evlo,evla,evdp,mag,stanum])
    elif ch == 'PandZ':
        # evttime, evlo, evla, evdp, mag, trnum_common
        evt_lst_P = os.listdir(ch_path['P'])
        evt_lst_Z = os.listdir(ch_path['Z'])
        evt_lst = evt_lst_P[:]
        for evt in evt_lst_P:
            if not evt in evt_lst_Z:
                evt_lst.remove(evt)
        for evt in evt_lst:
            path_evt_P = os.path.join(ch_path['P'],evt)
            path_evt_Z = os.path.join(ch_path['Z'],evt)
            sac_lst_P = [os.path.basename(sac).replace('.P','') for sac in glob.glob(os.path.join(path_evt_P,evt+'*.P'))]
            sac_lst_Z = [os.path.basename(sac).replace('.Z','') for sac in glob.glob(os.path.join(path_evt_Z,evt+'*.Z'))]
            sac_lst = sac_lst_P[:]
            for sac in sac_lst_P:
                if sac not in sac_lst_Z:
                    sac_lst.remove(sac)
            
            if len(sac_lst_P) == 0: os.rmdir(path_evt_P)
            elif len(sac_lst_Z) == 0: os.rmdir(path_evt_Z)
            elif len(sac_lst) == 0: print('No common traces for %s'%evt)
            else:    
                stats = obspy.read(os.path.join(path_evt_P,sac_lst_P[0]+'.P'))[0].stats.sac
                evlo, evla, evdp, mag = stats.evlo, stats.evla, stats.evdp, stats.mag
                stanum = len(sac_lst)
                evt_info.append([evt,evlo,evla,evdp,mag,stanum])

    elif ch == 'PorZ':
        # evttime, evlo, evla, evdp, mag, trnum_P, trnum_Z
        evt_lst_P = os.listdir(ch_path['P'])
        evt_lst_Z = os.listdir(ch_path['Z'])
        evt_lst = np.unique(evt_lst_P+evt_lst_Z).tolist()
        for evt in evt_lst:
            path_evt_P = os.path.join(ch_path['P'],evt)
            path_evt_Z = os.path.join(ch_path['Z'],evt)
            sac_lst_P = glob.glob(os.path.join(path_evt_P,evt+'*.P'))
            sac_lst_Z = glob.glob(os.path.join(path_evt_Z,evt+'*.Z'))
            if len(sac_lst_P) != 0 or len(sac_lst_Z) != 0:
                try:
                    stats = obspy.read(sac_lst_P[0])[0].stats.sac
                except:
                    stats = obspy.read(sac_lst_Z[0])[0].stats.sac
                evlo, evla, evdp, mag = stats.evlo, stats.evla, stats.evdp, stats.mag
                stanum_P, stanum_Z = len(sac_lst_P), len(sac_lst_Z)
                evt_info.append([evt,evlo,evla,evdp,mag,stanum_P,stanum_Z])
    
    else:
        print('wrong channel!')
        return

    with open(os.path.join(evtlog_path,'evt_%s.log'%ch),'w+') as f:
        for evtline in evt_info:
            f.write(' '.join([str(tmp) for tmp in evtline])+'\n')
        f.close()

# read config
config = configparser.ConfigParser()
parafile = 'paraP.cfg'
if os.path.isfile(parafile): config.read(parafile)
else: Usage()
ch_path = {}
out_path = config.get('path','out_path')
Z_path = config.get('path','Z_dir')
ch_path['Z'] = os.path.join(out_path,Z_path)
P_path = config.get('path','P_dir')
ch_path['P'] = os.path.join(out_path,P_path)

for ch in ['P','Z','PandZ','PorZ']:
    evt_catalog(ch,ch_path,evtlog_path)
