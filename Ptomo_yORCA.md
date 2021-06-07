# Ptomo_yORCA

## I.Differential travel time measurement

### 1. Data collection (./, no need to redo that for yORCA)

Collect segments including *P*-wave signal

(1) config file: ***paraP.cfg*** the default time interval is [-70, 30]s (relative to taup *P* arrival)

(2) data collecting script: ***collect_P.py***

Usage: *python collect_P.py -Sstation -Cchannel[all|P|Z] para.cfg*

parallel version: ***python parallel_collect.py***

output: *P*-wave signal sac file for pressure and vertical channel respectively in ***Pdata/*** like *Pdata/pressure/20190217112917/20190217112917_WW04.P*

(3) Event catalog (necessary before doing MCCC)

***python evtlog.py***

output: evtlog for Z, P, ZandP, ZorP in ***evtlog/***

***note: please don't change current event catalog for yORCA measurement, that means you shouldn't run python evtlog.py***

(4) azimuthal distribution plot (optional)

plot all candidate events

***python event_map.py -Ssta[all|sta] -C[P|Z|logfile] para.cfg***

output: evtmap.png 

only plot events used in tomography

***python event_map_log.py***

output: evtmap_log.png 



### 2. GUI program for MCCC (cc/)

GUI Program: ***cc_gui.py***

the config (head) file and channel can be specified in 

*line 48 def get_para(head='../paraP.cfg',ccch='P'):*

​	*\# ccch: P, Z, PandZ, PorZ*

some parameters related to interactive MCCC can be adapted in the main() function

run ***python cc_gui.py***

Input: sac file in ***../Pdata/***

output:

**$dT$ result format:** evtopts.evttime, rayp_rad, gcarc, baz, sta, nwk, orid, evla, evlo, evdp, dt, dcstd, dvcstd, acor, rsstd (std from resampling), cf (central frequency), ch

**file:** *cclog/{ch}\_mag\_{magmin}\_{magmax}/{ch}\_mag\_{magmin}\_{magmax}.dat*

**window/frequency record format:** left, right, freqmin, freqmax

**file:** *cclog/{ch}\_mag\_{magmin}\_{magmax}/{ch}\_mag\_{magmin}\_{magmax}\_winfreq.dat*



(1) Typical procedure for doing MCCC 

adjust frequency band $\rightarrow$ ampcull  $\rightarrow$ deleteBadTrs $\rightarrow$ tap '$\uparrow$' or '$\downarrow$' on the keyboard to adjust waveform amplitude displayed  $\rightarrow$ click bad traces to reject them $\rightarrow$ adjust time interval for P arrival $\rightarrow$ click 'reMCCC' to align and sort the segment $\rightarrow$ further culling and reMCCC  $\rightarrow$ tap 's' in the keyboard to save the measurements into cache (-> click 'finish' to save the measurements into the measurement dat file, you also can do this after doing MCCC for multiple events) $\rightarrow$ click 'Next(N)' to process the next event

(2) Changing measurements based on existing measurement dat file

everytime we run cc_gui.py, the script would automatically load the existing measurements in **cclog/{ch}\_mag\_{magmin}\_{magmax}/**, if you redo MCCC for a event with different parameters and want to save the new result, just tap 's' in the keyboard after doing MCCC to save the change as well as the new MCCC figure, don't forget to click 'finish' to save those change to a dat file.

***note: for the yORCA data, I suggest to divide the events into Mw5.5-6 and Mw6-10 for vertical records, and only use Mw6-10 for pressure records, you can change the magnitude range in the main() function of cc_gui.py***



### 3. measurement combination (cc/)

The basic idea is to combine **Z_mag_6.0_10.0.dat** and **P_mag_6.0_10.0.dat**, then concatenate the combined data file and **Z_mag_5.5_6.0.dat**

(1) direct combination

***python combinePZ.py***

you need to specify the parameter **dominant**

**dominant = 'Z' or 'P':*** only retain vertical or pressure measurement for the common event

**dominant = 'num':** only retain the channel with more traces for the common event

**dominant = 'both':** retain all measurements from both channel but with 2 different event index

output: **PZ_direct_{dominant}.dat**

(2) mixed combination

***python ZPmccc.py***

No need to specify any parameters

output: **PZmix.dat**



### 4. dT distribution plot

you need to specify the dT dat_file in dtresult.py

and then run ***python dtresult.py***

output: dt distribution figure (.png)



## II. Tomography (tomo/BWTOMOG_atten_Vp_Vs_joint/)

config file: ***yORCA/PARMS.m***

data file: *yORCA/data/{Z_mag_5.5_10.0.dat} {P_mag_6.0_10.0.dat} {PZ_direct_both.dat} {PZ_direct_num.dat} {PZmix.dat}*

sta file: *yORCA/data/stations_PZ.dat*

specify data file and regularization parameters in ***PARMS.m*** and run ***RUN_INVERSION_ORCA.m*** for 3D tomographic inversion

the model files were saved in ***model/***

you can load them and use ***modelplot.m*** to plot their tomographic figures

### 1. Ltest

firstly run ***RUN_LTEST.m***

you should specify data file in ***RUN_LTEST.m***

output: Ltest_flatness/Ltest_dT_*.mat

then run ***plot_Lcurves.m***

you should specify the figure saving path, the smooth type (flatness/smoothness), the F value (resid2rough) in ***plot_Lcurves.m***

e.g. ofile = 'Ltest_dT_PZmix.mat'; %0.6 is the desired F value

output: Lcurve plot and contour plot of penalty function saved in specified address



### 2. checkerboard test

use the function syn_test_ORCA(syn_datafile,damp,smooth) defined in ***synth_test_ORCA.m***



### 3. squeezing test

(1) one-side damped squeezing

use the function squeezing_test(syn_datafile,direction,is2D) defined in ***squeezing_test.m***

direction: 'up' (downside-damped) or 'down' (upside-damped)

is2D: 0 (3D inversion) or 1 (2D inversion)

output: ***sqztest/sqztest_\*.mat*** (sqztest numerical result), tomographic figures

then use ***sqzplot.m*** to plot the variation curve of variance reduction and m2 norm



(2) 2-side damped squeezing

use the function sqztest_2side(syn_datafile,is2D) defined in ***sqztest_2side.m***

output: ***sqztest_2side/sqz2side_\*.mat***  (2 side damped sqztest numerical result)



(3) move scan with a depth window containing a certain number of continuous layer

use ***sqztest_2side_movescan.m***, you will need to specify the number of layers in the depth window

output: the plot showing the variance reduction and m2 norm



### 4. 2D inversion

(1) find optimal strike direction of the profile

run ***RUN_ROTATE_TEST.m***, please specify the data file used in ***RUN_ROTATE_TEST.m***

output: ***rotate_test_\*.mat***

then run ***rotatetestplot.m***

output: plots for data misfit vs. strike direction



(2) 2D tomographic inversion

in yORCA/PARMS.m

on line 80, change par.dampx from 1 to 1000

on line 117, change mstruct.origin from [par.origin 0] to [par.origin -15]

then run ***RUN_INVERSION_ORCA.m*** for 2D inversion