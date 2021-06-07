#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 08:16:59 2020

@author: lun
"""
import numpy as np
from scipy.io import loadmat
from scipy.interpolate import griddata
import mayavi
from mayavi import mlab

def mesh4D(x,y,z,v,extent=[]):
    x,y,z,v = x.copy(),y.copy(),z.copy(),v.copy()
    if len(extent) == 6:
        ind = reduce(np.intersect1d,(np.where(extent[0]<=x)[0],np.where(x<=extent[1])[0],\
                np.where(extent[2]<=y)[0],np.where(y<=extent[3])[0],\
                np.where(extent[4]<=z)[0],np.where(z<=extent[5])[0]))
        x,y,z,v = x[ind],y[ind],z[ind],v[ind]
    xuni,yuni,zuni = np.array(sorted(np.unique(x))),np.array(sorted(np.unique(y))),np.array(sorted(np.unique(z)))
    xdim,ydim,zdim = len(xuni),len(yuni),len(zuni)
    dx = float(np.diff(xuni)[0])
    dy = float(np.diff(yuni)[0])
    dz = float(np.diff(zuni)[0])
    xmin,ymin,zmin = min(xuni),min(yuni),min(zuni)
    xmax,ymax,zmax = max(xuni),max(yuni),max(zuni)
    vgrid = np.zeros([xdim,ydim,zdim])
    for i in range(v.shape[0]):
        i1,i2,i3 = np.where(xuni==x[i])[0][0],np.where(yuni==y[i])[0][0],np.where(zuni==z[i])[0][0]
        vgrid[i1,i2,i3] = v[i]
    #xmesh,ymesh,zmesh = np.meshgrid(x,y,z)
    xmesh,ymesh,zmesh = np.mgrid[xmin:xmax+dx:dx,ymin:ymax+dy:dy,zmin:zmax+dz:dz]
    #xmesh,ymesh,zmesh = np.meshgrid(xuni,yuni,zuni)
    return xmesh,ymesh,zmesh,vgrid

G_0 = loadmat('G_0.mat')['G_0'].toarray()
par = loadmat('par.mat')['par']
mx = np.array(par['mx'][0,0]).reshape(1,-1)[0]
my = np.array(par['my'][0,0]).reshape(1,-1)[0]
mz = np.array(par['mz'][0,0]).reshape(1,-1)[0]

xmax, ymax, zmax = map(np.max,[mx,my,mz])
xmin, ymin, zmin = map(np.min,[mx,my,mz])
extent = [xmin,xmax,ymin,ymax,zmin,zmax]


xg, yg, zg = np.mgrid[xmin:xmax:20,ymin:ymax:20,zmin:zmax:20]
xg = xg.reshape(1,-1)[0]
yg = yg.reshape(1,-1)[0]
zg = zg.reshape(1,-1)[0]
K = griddata(np.vstack((mx,my,mz)).T,G_0[4],np.vstack((xg,yg,zg)).T,method='linear')
# 0.1,0.3,0.5,1,2,5,10

xmesh,ymesh,zmesh,vgrid = mesh4D(xg,yg,zg,K)

m = mlab.contour3d(xmesh,ymesh,zmesh,vgrid,contours=100)
#mlab.points3d(mx,my,mz,G_0[0])
mlab.axes(xlabel='x[km]', ylabel='y[km]', zlabel='Depth[km]', nb_labels=3)
mlab.colorbar(m,orientation='vertical',nb_labels=5,title='sensitivity')
mlab.show()