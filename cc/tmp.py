import numpy as np
import matplotlib.pylab as plt
Zdat = 'Z_mag_5.5_10.0.dat'
Pdat = 'P_mag_6.0_10.0.dat'
PZmixdat = 'PZmix.dat'
evt_time = '20180916211148'
Zevtid = np.loadtxt(Zdat,usecols=0,dtype='str')
Pevtid = np.loadtxt(Pdat,usecols=0,dtype='str')
PZevtid = np.loadtxt(PZmixdat,usecols=0,dtype='str')
Zdcor = np.loadtxt(Zdat,usecols=10)
Pdcor = np.loadtxt(Pdat,usecols=10)
PZdcor = np.loadtxt(PZmixdat,usecols=10)
Zstalst = np.loadtxt(Zdat,usecols=4,dtype='str')
Pstalst = np.loadtxt(Pdat,usecols=4,dtype='str')
PZstalst = np.loadtxt(PZmixdat,usecols=4,dtype='str')

allsta = PZstalst[PZevtid==evt_time].tolist()
Zsta = Zstalst[Zevtid==evt_time].tolist()
Psta = Pstalst[Pevtid==evt_time].tolist()
staid = {}
allsta.remove('EC02')
for i in range(len(allsta)):
    staid[allsta[i]] = i
staid['EC02'] = i+1
allsta += ['EC02']

fig = plt.figure(figsize=(14,5))
ax = fig.add_subplot(111)
i = -0.15
k = 0
label = np.zeros(3)
for evtid, dcorlst, stalst, ch, color in zip([Zevtid,Pevtid,PZevtid],[Zdcor,Pdcor,PZdcor],[Zstalst,Pstalst,PZstalst],['Z','P','PZmix'], ['b','r','g']):
    substalst = stalst[evtid==evt_time]
    subdcor = dcorlst[evtid==evt_time]
    for sta, dcor in zip(substalst,subdcor):
        if label[k] == 0:
            ax.vlines(staid[sta]+i,0,dcor,color,linewidth=5,label=ch)
            label[k] = 1
        else:
            ax.vlines(staid[sta]+i,0,dcor,color,linewidth=5)
    i += 0.15
    k += 1
ax.set(xticks=np.arange(0,len(allsta)),xticklabels=allsta,xlim=[-0.5,len(allsta)-0.5],ylim=[-0.65,0.65],ylabel='dcor[s]')
ax.grid()
ax.legend(loc='upper right')
plt.suptitle(evt_time)
plt.savefig(evt_time+'.png')
plt.show()
