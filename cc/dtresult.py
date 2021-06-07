import numpy as np
import matplotlib.pylab as plt
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata
from scipy.io import loadmat
import os, glob
#import distaz

def spoke(ax,az,origin=[0,0],length=1,color='k',linewidth=2):
    dy = length*np.cos(np.deg2rad(az))
    dx = length*np.sin(np.deg2rad(az))
    ax.plot([origin[0],origin[0]+dx],[origin[1],origin[1]+dy],color=color,linewidth=linewidth,zorder=20)

#DCOR_path = '../NOISETC_CI/DATA/NOISETC/DCOR/DCOR(1.0to10.0s)'

#DCOR_path = '../NOISETC_CI/DATA/NOISETC/DCOR/DCOR(1.0to10.0s)'
#mseed_path = '/home/lun/scratch-lun/OBSdata/yORCA_data/Mseed'
sta_latlon = '../../background/latlon_all.dat'
GMT_dir = 'GMT_dir'
sta_info = {}
sta_dt = {}
#dat_file = 'Z_mag_5.5_10.0.dat'
#dat_file = 'P_mag_6.0_10.0.dat'
#dat_file = 'PZmix.dat'
dat_file = 'PZ_direct_num.dat'

with open(sta_latlon,'r') as f:
    for line in f.readlines():
        line = line.replace('\n','').split()
        sta0, stlo0, stla0, stdp0 = line[0], float(line[1]), float(line[2]), float(line[3])
        sta_info[sta0] = [stlo0, stla0, stdp0] #lon, lat, dep
        sta_dt[sta0] = []
    f.close()

stacol = np.loadtxt(dat_file,usecols=4,dtype='str')
bazcol = np.loadtxt(dat_file,usecols=3)
dtcol = np.loadtxt(dat_file,usecols=10)
stdcol = np.loadtxt(dat_file,usecols=11)

stalst = np.unique(stacol).tolist()

cNorm = colors.Normalize(vmin=-0.5, vmax=0.5)
scalarMap = plt.cm.ScalarMappable(norm=cNorm, cmap=plt.cm.rainbow)
fig = plt.figure(figsize=(6,5))
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
figxb = 0.1
figyb = 0.1
figlen = 0.75
figwid = 0.75
annote_size = 9
#plt.axis("equal")
ax = fig.add_axes([figxb,figyb,figlen,figwid])
ax.axis('equal')
for item in (ax.get_xticklabels()+ax.get_yticklabels()): item.set_fontsize(annote_size)

#R = '-137/-129/-10/-3'
#xb, xe, yb, ye = -136, -130, -10, -3
xb, xe, yb, ye = -140, -126, -11.5, -1
dx, dy = 0.05, 0.05
#os.system('gmt grdfft %s/gravity1.grd -Fr5/2.5/0.5/0.3 -G%s/tmp1.grd'%(GMT_dir,GMT_dir))
#os.system('gmt grd2xyz %s/tmp1.grd > %s/SPsample.txt'%(GMT_dir,GMT_dir))
txt = np.loadtxt('%s/SPsample.txt'%GMT_dir)

#xb, xe, yb, ye = np.min(txt[:,0]), np.max(txt[:,0]), np.min(txt[:,1]), np.max(txt[:,1])
dim1 = int((ye-yb)/dy)+1
dim2 = int((xe-xb)/dx)+1
lat = np.linspace(ye,yb,num=dim1,endpoint=True)
lon = np.linspace(xe,xb,num=dim2,endpoint=True)
A = np.zeros([dim1,dim2])

lon0, lat0 = np.linspace(xb,xe,100), np.linspace(yb,ye,100)
grid_lon, grid_lat = np.meshgrid(lon0,lat0)
grid_grav = griddata(txt[:,1::-1],txt[:,2],(grid_lat,grid_lon),method='nearest')
print(grid_grav)

#ax.contour(lon0,lat0,grid_grav,cmap=plt.cm.gist_gray,levels=np.linspace(-20,20,num=10,endpoint=True),extent=[xb,xe,yb,ye],linestyles='solid',linewidths=0.5)
#ax.set(xlim=[xb,xe],ylim=[yb+1,ye-1])
xticks, yticks = np.arange(-136,-129,1), np.arange(-9,-3,1)
xticklabels, yticklabels = [], []
for xtick in xticks: xticklabels.append('$%d^{\circ}$'%xtick)
for ytick in yticks: yticklabels.append('$%d^{\circ}$'%ytick)
ax.set(xlim=[-136,-130],ylim=[-9,-3.5],xticks=xticks,xticklabels=xticklabels,\
    yticks=yticks,yticklabels=yticklabels)
cb = ax.contourf(lon0,lat0,grid_grav,cmap=plt.cm.Greys,\
    levels=np.linspace(-20,20,num=11,endpoint=True))
c_ax1 = fig.add_axes([figxb,figyb+figwid+0.02,figlen,0.025])
fig.colorbar(cb,cax=c_ax1,\
    orientation='horizontal',label='$free-air gravity(mgal)$')
c_ax1.xaxis.set_ticks_position('top')
c_ax1.xaxis.set_label_position('top')
for item in c_ax1.get_xticklabels(): item.set_fontsize(annote_size)   


for sta,dt,baz,std in zip(stacol,dtcol,bazcol,stdcol):
    colorVal = scalarMap.to_rgba(dt)
    stlo, stla = sta_info[sta][:2]
    #spoke(ax,baz,origin=[stlo,stla],length=0.25,color=colorVal,linewidth=std*20)
    spoke(ax,baz,origin=[stlo,stla],length=0.25,color=colorVal)
    sta_dt[sta].append(dt)

for sta in stacol:
    #if sta == 'CC02':
    #    continue
    if len(sta_dt[sta]) == 0:
        continue
    stlo, stla, stdp = sta_info[sta]
    print(sta,np.min(sta_dt[sta]),np.max(sta_dt[sta]))
    avg_dt = np.mean(sta_dt[sta])
    colorVal = scalarMap.to_rgba(avg_dt)
    ax.scatter(stlo,stla,marker='o',s=70,c=colorVal,edgecolors='black',linewidths=1.5,zorder=30)

ax.grid(linestyle='--',linewidth=0.5)
ax.set(xlabel='$Longitude$',ylabel='$Latitude$')
print('Total P arrival: %d'%(len(stacol)))
scalarMap.set_array([])
c_ax2 = fig.add_axes([figxb+figlen+0.02,figyb,0.02,figwid])
fig.colorbar(scalarMap,cax=c_ax2,label=r'$\longleftarrow FAST$     $\Delta T(s)$     $SLOW \longrightarrow$')
for item in c_ax2.get_yticklabels(): item.set_fontsize(annote_size)
plt.savefig('dtfinal_PZ_direct_num%d.png'%len(stacol))
plt.show()
