import sys, getopt
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pylab as plt
import distaz, obspy, configparser, os, glob, geo, gc

sta = 'all'
logfile = 'tomo/BWTOMOG_atten_Vp_Vs_joint/yORCA/data/PZ_direct_num.dat'
stapos_lst = '../background/latlon_all.dat'

def dist_contour(lon0,lat0,gcarc):
    azi = np.linspace(0,360,num=100,endpoint=True)
    lonlat = np.zeros([len(azi),2])
    for i in range(len(azi)):
        lat, lon = geo.latlon_from(lat0,lon0,azi[i],np.array([gcarc]))
        lonlat[i,0], lonlat[i,1] = lon[0], lat[0]
    return lonlat

def count_evt_logfile(logfile):
    evt_info = []
    evtid = np.loadtxt(logfile,usecols=6)
    evla = np.loadtxt(logfile,usecols=7)
    evlo = np.loadtxt(logfile,usecols=8)
    evdp = np.loadtxt(logfile,usecols=9)
    evtid_uni = np.unique(evtid)
    for i in range(len(evtid_uni)):
        ind = np.where(evtid==evtid_uni[i])[0][0]
        evt_info.append([evlo[ind],evla[ind],evdp[ind],0,sum(evtid==evtid_uni[i])])
    evt_info = np.array(evt_info)
    return evt_info

def make_polar_bar(ax,clon,clat,lonlst,latlst,
                        edgecolor='white',
                        bottom=0,
                        bins=18,
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
    ytick = np.arange(0,max(freq),5)
    yticklabel = ['%d'%tmp for tmp in ytick]
    
    ax.set(xticks=xtick,xticklabels=xticklabel,\
        ylim=[0,max(freq)],yticks=ytick,yticklabels=yticklabel)

lonlat_evt = count_evt_logfile(logfile)

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
plt.savefig('evtmap_log.png')
plt.show()
