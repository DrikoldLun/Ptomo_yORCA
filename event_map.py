import sys, getopt
def Usage():
    print('python event_map.py -Ssta[all|sta] -C[P|Z|logfile] para.cfg')
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
        if not ch_use in ['P','Z']:
            print("Channel needs to be 'P' or 'Z'!")
            Usage()
    else:
        Usage()

from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pylab as plt
import distaz, obspy, configparser, os, glob, geo, gc

config = configparser.ConfigParser()
parafile = sys.argv[-1]
if os.path.isfile(parafile): config.read(parafile)
else: Usage()
ch_path = {}
out_path = config.get('path','out_path')
Z_path = config.get('path','Z_dir')
ch_path['Z'] = os.path.join(out_path,Z_path)
P_path = config.get('path','P_dir')
ch_path['P'] = os.path.join(out_path,P_path)
stapos_lst = config.get('path','stapos_lst')

def dist_contour(lon0,lat0,gcarc):
    azi = np.linspace(0,360,num=100,endpoint=True)
    lonlat = np.zeros([len(azi),2])
    for i in range(len(azi)):
        lat, lon = geo.latlon_from(lat0,lon0,azi[i],np.array([gcarc]))
        lonlat[i,0], lonlat[i,1] = lon[0], lat[0]
    return lonlat

def count_evt(sta,ch):
    evt_lst = os.listdir(ch_path[ch])
    if sta == 'all':
        evt_info = [] # evlo, evla, evdp, mag, stanum
        for evt in evt_lst:
            path_evt = os.path.join(ch_path[ch],evt)
            sac_lst = glob.glob(os.path.join(path_evt,evt+'*.'+ch))
            if len(sac_lst) == 0: os.rmdir(path_evt)
            else:
                stats = obspy.read(sac_lst[0])[0].stats.sac
                evlo, evla, evdp, mag = stats.evlo, stats.evla, stats.evdp, stats.mag
                stanum = len(sac_lst)
                evt_info.append([evlo,evla,evdp,mag,stanum])
                del stats
                gc.collect()
    else:
        evt_info = [] # evlo, evla, evdp, mag
        for evt in evt_lst:
            path_evt = os.path.join(ch_path[ch],evt)
            sac = os.path.join(path_evt,evt+'_'+sta+'.'+ch)
            if not os.path.exists(sac): pass
            else:
                stats = obspy.read(sac)[0].stats.sac
                evlo, evla, evdp, mag = stats.evlo, stats.evla, stats.evdp, stats.mag
                evt_info.append([evlo,evla,evdp,mag])
                del stats
                gc.collect()
    evt_info = np.array(evt_info)
    return evt_info

def make_polar_bar(ax,clon,clat,lonlst,latlst,
                        edgecolor='white',
                        bottom=0,
                        bins=36,
                        opening=1,
                        ticks=None,
                        figsize=(8, 8),
                        isplot=False,
                        **kwargs):

    bazi = np.zeros(len(lonlst))
    i = 0
    for lon0,lat0 in zip(lonlst,latlst):
        bazi[i] = distaz.DistAz(clat,clon,lat0,lon0).getBaz()
        i += 1
    #print(bazi)
    theta  = np.linspace(0.0,360,bins,endpoint=True)
    freq, theta = np.histogram(bazi, bins=theta)
    theta = np.deg2rad((theta[:-1]+theta[1:])/2)
    width = (2 * np.pi) / bins * opening

    bars = ax.bar(
        theta,
        freq,
        width=width,
        bottom=bottom,
        edgecolor='white',
        facecolor='black',
        **kwargs)
    '''
    for f,bar in zip(freq,bars):
        bar.set_facecolor(plt.cm.viridis(float(f)/len(array_data)))
        bar.set_alpha(0.5)
    '''
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    xtick = np.deg2rad(np.arange(0,360,45))   
    xticklabel = ['N', '$45^{\circ}$', 'E', '$135^{\circ}$', 'S', '$225^{\circ}$', 'W', '$315^{\circ}$']
    ytick = np.arange(0,max(freq),50)
    yticklabel = ['%d'%tmp for tmp in ytick]
    
    ax.set(xticks=xtick,xticklabels=xticklabel,\
        ylim=[0,max(freq)],yticks=ytick,yticklabels=yticklabel)

lonlat_evt = count_evt(sta,ch_use)

if sta == 'all':
    lonlat_sta = np.loadtxt(stapos_lst,usecols=(1,2))
    [lon0,lat0] = lonlat_sta.mean(axis=0)
else:
    with open(stapos_lst,'r') as f:
        for line in f.readlines():
            line = line.replace('\n','').split()
            if line[0] == sta:
                lon0, lat0 = [float(tmp) for tmp in line[1:3]]
                break
        f.close()
print(lon0,lat0)

contour30 = dist_contour(lon0,lat0,30)
contour95 = dist_contour(lon0,lat0,96)

fig = plt.figure(figsize=(10,5))
ax1 = fig.add_axes([0.1,0.1,0.35,0.8])

ax1.set_title('Event map',pad=25)
m = Basemap(projection='aeqd',lon_0=lon0,lat_0=lat0,resolution='l')
#m.etopo()
m.shadedrelief()
x, y = m(lonlat_evt[:,0],lonlat_evt[:,1])
x0, y0 = m(lon0,lat0)
m.plot(x0,y0,'y^',markersize=8)
m.plot(contour30[:,0],contour30[:,1],'k--',latlon=True,linewidth=1)
m.plot(contour95[:,0],contour95[:,1],'k--',latlon=True,linewidth=1)
if sta == 'all':
    cb = m.scatter(x,y,marker='o',s=30,c=lonlat_evt[:,4],\
            cmap=plt.cm.rainbow,edgecolors='none',zorder=100)
    c_ax = fig.add_axes([0.46,0.2,0.01,0.6])
    fig.colorbar(cb,cax=c_ax,label='station number')
else:
    m.scatter(x,y,marker='o',s=30,color='none',edgecolors='red',zorder=100)

ax2 = fig.add_axes([0.55,0.1,0.35,0.8],polar='True')
make_polar_bar(ax2,lon0,lat0,lonlat_evt[:,0],lonlat_evt[:,1])
ax2.set_title('Azimuth distribution (num=%d)'%np.size(lonlat_evt,0),pad=10)
plt.savefig('evtmap_%s.png'%ch_use)
plt.show()
