import sys, getopt
def Usage():
    print('python collect_P.py -Sstation -C[all|P|Z] para.cfg')
    sys.exit(1)
try:
    opts, args = getopt.getopt(sys.argv[1:],'S:C:')
except:
    Usage()
if opts == []:
    Usage()

ch_use = 'all'
for op, value in opts:
    if op == '-S':
        sta = value
    elif op == '-C':
        ch_use = value
        if not ch_use in ['all','P','Z']:
            print("Channel needs to be 'all', 'P' or 'Z'!")
            Usage()
    else:
        Usage()

import obspy,os,configparser,glob
import numpy as np
from datetime import date
from obspy.taup import TauPyModel
import distaz, geo

config = configparser.ConfigParser()
parafile = sys.argv[-1]
if os.path.isfile(parafile): config.read(parafile)
else: Usage()

data_path = config.get('path','data_path')
out_path = config.get('path','out_path')
Z_path = config.get('path','Z_dir')
Z_path = os.path.join(out_path,Z_path)
P_path = config.get('path','P_dir')
P_path = os.path.join(out_path,P_path)
image_path = config.get('path','image_dir')
image_path = os.path.join(out_path,image_path)
resp_path = config.get('path','resp_path')
evt_lst = config.get('path','evt_lst')
stapos_lst = config.get('path','stapos_lst')
culling_path = config.get('path','culling_path')
for tmp_path in [Z_path,P_path,image_path]:
    if not os.path.exists(tmp_path): os.makedirs(tmp_path)
time_beforetauP = config.getfloat('para','time_beforetauP') 
time_aftertauP = config.getfloat('para','time_aftertauP')
dismin = config.getfloat('para','dismin')
dismax = config.getfloat('para','dismax')
magmin = config.getfloat('para','magmin')
samprate = config.getfloat('para','samprate')
#noisegate = config.getfloat('para','noisegate')

def preprocess(seis):
    seistmp = seis.copy()
    seistmp.detrend(type="demean")
    seistmp.detrend(type="linear")
    seistmp.taper(max_percentage=0.05)
    sta = seistmp.stats.station
    ch = seistmp.stats.channel
    respfile = os.path.join(resp_path,ch,'RESP.XX.'+sta+'..'+ch)
    pre_filt = (0.005, 0.01, 20.0, 25.0)
    seedresp = {'filename':respfile,
                'units':'VEL'}
    seistmp.simulate(pre_filt=pre_filt,seedresp=seedresp)
    seistmp.detrend(type='constant')
    seistmp.detrend(type='linear')
    return seistmp

with open(evt_lst,'r') as f:
    event_lst = f.readlines()
    f.close()

stlo, stla, stdp = {}, {}, {}
with open(stapos_lst,'r') as f:
    for line in f.readlines():
        line = line.replace('\n','').split()
        sta0, stlo0, stla0, stdp0 = line[0],float(line[1]), float(line[2]), float(line[3])
        stlo[sta0] = stlo0
        stla[sta0] = stla0
        stdp[sta0] = stdp0

model = TauPyModel(model="iasp91")

def collect_sac(sta):
    # data culling
    culling_hour = []
    for CH in ['CH2','CH3']:
        with open(os.path.join(culling_path,sta+'_'+CH+'_culling.dat'),'r') as f:
            for line in f.readlines():
                line = line.replace('\n','').split()
                culling_hour.append([CH]+[int(tmp) for tmp in line[2:6]])
            f.close()

    evnum_Z, evnum_P = 0, 0
    for event in event_lst:
        event = event.replace('\n','').split()
        #hour_index = float(event[0])
        evlo = float(event[3])
        evla = float(event[4])
        evdp = float(event[5])
        mag = float(event[1])
        tb = obspy.UTCDateTime(event[2])
        y, m, d = tb.year, tb.month, tb.day
        d_year = (date(y,m,d)-date(y,1,1)).days+1
        h_event = tb.hour
        min_event = tb.minute
        sec_event = tb.second
        
        evt_time_txt = '%04d%02d%02d%02d%02d%02d'%(y,m,d,h_event,min_event,sec_event)
        # reject events based on magnitude
        if mag < magmin:
            print(sta+': '+evt_time_txt+' mag %.1f < threshold %.1f'%(mag,magmin))
            continue
        result = distaz.DistAz(stla[sta],stlo[sta],evla,evlo)
        gcarc = result.getDelta()
        az = result.getAz()
        baz = result.getBaz()
        dist = result.degreesToKilometers(gcarc)
        # reject events based on distance
        if not (dismin <= gcarc <= dismax):
            print(sta+': '+evt_time_txt+' gcarc %.1f is outside of bound (%.1f,%.1f)'%(gcarc,dismin,dismax))
            continue
        try:
            arrivals = model.get_travel_times(source_depth_in_km=evdp,distance_in_degree=gcarc,phase_list=['P'])
            t_P = arrivals[0].time
            rayp = geo.srad2skm(arrivals[0].ray_param)
        except:
            print(sta+': '+evt_time_txt+' no P arrival')
            continue
        
        tmp_path = os.path.join(data_path,sta,sta+'.*.'+str(y)+'.'+'%03d'%(d_year)+'*.msd')
        if len(glob.glob(tmp_path)) == 0:
            print(sta+': '+evt_time_txt+' no msd file')
            continue

        new_ch = ['Z','P']
        ch_path = [Z_path,P_path]
        if ch_use == 'all':
            b_ch, e_ch = 0, 2
        elif ch_use == 'Z':
            b_ch, e_ch = 0, 1
        elif ch_use == 'P':
            b_ch, e_ch = 1, 2
        for i in range(b_ch,e_ch):
            CH = 'CH'+str(i+2)
            if [CH,y,m,(tb+t_P).hour] in culling_hour:
                print(sta+': '+evt_time_txt+' '+new_ch[i]+' in culling hour')
                continue
            msd = os.path.join(data_path,sta,'.'.join([sta,CH,str(y),'%03d'%d_year])+'*.msd')
            if len(glob.glob(msd)) == 0: continue
            st = obspy.read(msd)[0].copy()
            st.data = st.data.astype('float64')
            fs = int(st.stats.sampling_rate)
            diff = tb-st.stats.starttime
            #if diff < 0:
            #    print(sta+': '+evt_time_txt+' '+new_ch[i]+' discontinuity between 2 days')
            #    continue
            nb = int((diff+t_P-time_beforetauP)*fs)
            ne = int((diff+t_P+time_aftertauP)*fs+1)
            if nb < 0 or ne > len(st.data):
                print(sta+': '+evt_time_txt+' '+new_ch[i]+' discontinuity between 2 days')
                continue
            st.data = st.data[nb:ne]
            st.stats.starttime = tb+t_P-time_beforetauP
            try:
                st = preprocess(st)
            except:
                print(sta + ": can't preprocess "+msd)
                continue
            evt_ch_path = os.path.join(ch_path[i],evt_time_txt)
            if not os.path.exists(evt_ch_path): os.mkdir(evt_ch_path)
            sac = os.path.join(evt_ch_path,evt_time_txt+'_%s.%1s'%(sta,new_ch[i]))
            st.stats.channel = new_ch[i]
            st.write(sac,'SAC')
            try:
                st = obspy.read(sac)[0]
                st.stats.sac.evdp = evdp
                st.stats.sac.evla = evla
                st.stats.sac.evlo = evlo
                st.stats.sac.stla = stla[sta]
                st.stats.sac.stlo = stlo[sta]
                st.stats.sac.stdp = stdp[sta]
                st.stats.sac.mag = mag
                st.stats.sac.gcarc = gcarc
                st.stats.sac.dist = dist
                st.stats.sac.az = az
                st.stats.sac.baz = baz
                st.stats.sac.user0 = tb
                st.stats.sac.user1 = time_beforetauP
                st.stats.sac.user2 = time_aftertauP
                st.stats.sac.user3 = t_P
                st.stats.sac.user4 = rayp
                st.write(sac,'SAC')
            except:
                print("can't read "+sac)
                os.system('rm '+sac)
                continue
            print(sac)
            if i == 0: evnum_Z += 1
            elif i == 1: evnum_P += 1
    print('%s: %d Z, %d P found'%(sta,evnum_Z,evnum_P))

collect_sac(sta)
