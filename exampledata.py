path = '../RF/RF_data/event_data'
evtplace = 'S. ALASKA'
evt = '20181130172929'
#evt = '20181118202546'
#evt = '20181102110115'
ch = 'Z'
win = [0,2500]
freq = [0.01,0.12]

import obspy,os,glob
import numpy as np
import matplotlib.pylab as plt
import distaz
from obspy.taup import TauPyModel
model = TauPyModel(model="iasp91")

def preprocess(tr,be=[0,3600],freq=[0.01,0.12],win=[0,2500]):
    npts = tr.stats.npts
    time = np.linspace(be[0],be[1],num=npts,endpoint=True)
    tr.filter('bandpass',freqmin=freq[0],freqmax=freq[1],zerophase=True)
    tr.data = tr.data[np.where((time>=win[0])&(time<=win[1]))]
    tr.taper(max_percentage=0.05,type='hann')
    return tr.data

dist_lst = []
sta_lst = []
data = []

for sta in os.listdir(path):
    if sta in ['CC09','CC01']: continue
    sac = os.path.join(path,sta,'%s_%s.%s'%(evt[:12],sta,ch))
    if not os.path.isfile(sac): continue
    st = obspy.read(sac)[0]
    stats = st.stats.sac
    result = distaz.DistAz(stats.stla,stats.stlo,stats.evla,stats.evlo)
    dist_lst.append(result.getDelta())
    data.append(preprocess(st,be=[0,3600],freq=freq,win=win))
    sta_lst.append(sta)

dist_axis = np.linspace(min(dist_lst)-1,max(dist_lst)+1,100,endpoint=True)
phase_list = ['P','pP','PP','S','sS','SS','SSS']
phase_time = {}
phase_dist_list = {}
for phase in phase_list:
    phase_time[phase] = []
    phase_dist_list[phase] = []
    for dist in dist_axis:
        try:
            arrivals = model.get_travel_times(source_depth_in_km=stats.evdp,distance_in_degree=dist,phase_list=[phase])
            phase_time[phase].append(arrivals[0].time)
            phase_dist_list[phase].append(dist)
        except:
            continue

#arrivals = model.get_travel_times(source_depth_in_km=stats.evdp,distance_in_degree=np.mean(dist_lst),phase_list=['P','pP','PP','S','sS','SS','SSS']) 

time = np.linspace(win[0],win[1],num=len(data[0]),endpoint=True)
fig = plt.figure(figsize=(9,12))
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
ax = fig.add_subplot(111)
ax2 = ax.twinx()
for i in range(len(sta_lst)):
    ax.plot(time,data[i]/np.abs(data[i]).max()/2.+dist_lst[i],'k-',linewidth=1)
ax.set(xlim=win,xlabel='Time since earthquake (s)',ylim=[min(dist_lst)-1,max(dist_lst)+1],ylabel='Distance from Earth (deg)',title=evtplace+' Mw %.1f %4s-%2s-%2s BHZ (%.2f-%.2f Hz)'%(stats.mag,evt[:4],evt[4:6],evt[6:8],freq[0],freq[1]))
ax2.set(xlim=win,ylim=[min(dist_lst)-1,max(dist_lst)+1],yticks=dist_lst,yticklabels=sta_lst)
ax.invert_yaxis()

i = 0
for phase in phase_list:
    ax.plot(phase_time[phase],phase_dist_list[phase],color='blue',linestyle='--',linewidth=0.5)
    ax.text(phase_time[phase][0]-10,min(dist_lst)-1+(0.015*(-1)**(i-1)+0.05)*(max(dist_lst)-min(dist_lst)+2),phase,rotation=-90,fontstyle='italic',fontweight='bold',color='blue')
    i += 1
'''
for i in range(7):
    try:
        ax.vlines(arrivals[i].time,min(dist_lst)-1,max(dist_lst)+1,color='blue',linestyle='--',linewidth=0.5)
        ax.text(arrivals[i].time-10,max(dist_lst)+1-(0.03*(-1)**(i-1)+0.05)*(max(dist_lst)-min(dist_lst)+2),arrivals[i].name,rotation=-90,fontstyle='italic',fontweight='bold',color='blue')
    except:
        continue
'''
plt.savefig('yORCAexample.png')
plt.show()
