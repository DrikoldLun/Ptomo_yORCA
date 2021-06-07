#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 01:02:16 2020

@author: lun
"""
import numpy as np
import matplotlib.pylab as plt

datlst = ['Z_mag_5.5_10.0.dat','P_mag_6.0_10.0.dat','PZ_direct_num.dat','PZmix.dat']
title = ['Z(Mw>=5.5)','P(Mw>=6.0)','PZ_direct','PZ_mixed']

fig = plt.figure(figsize=(10,8))
for i in range(4):
    rsstd = np.loadtxt(datlst[i],usecols=-7) #-5 dvcstd -7 dcor -3 rsstd
    xhist, xedge = np.histogram(rsstd,bins=80)
    xprob = xhist/float(len(rsstd))
    xx = 0.5*(xedge[0:-1]+xedge[1:])
    ax = fig.add_subplot(2,2,i+1)
    width = xedge[1]-xedge[0]
    ax.bar(x=xx,height=xprob,width=width,facecolor='black',edgecolor='white')
    ax.plot(xx,xprob,'r-',label=title[i])
    ax.set(xlim=[-1,1],xlabel='dT[s]',ylabel='prob')
    ax.legend(loc='upper right')

plt.suptitle('dcor distribution')
plt.savefig('dcor.png')
plt.show()
