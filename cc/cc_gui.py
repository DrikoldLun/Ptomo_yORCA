#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 16:22:17 2020

@author: lun
"""
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 20:59:04 2020

@author: lun
"""

import matplotlib.pyplot as plt
import os, glob
import obspy
from obspy.signal.cross_correlation import correlate, xcorr_max
from scipy.signal import spectrogram
import numpy as np
from matplotlib.widgets import Button
from matplotlib.widgets import TextBox
import copy
import configparser
import cc
#import matlab.engine
#eng = matlab.engine.start_matlab()
config = configparser.RawConfigParser()
old12normal = ['CC02','CC04','CC05','CC07','CC08','WC03','WC04']
old12noise = ['CC06','EE04','CC11','WC05','WW03']
mid6noise = ['EE02']
mid6other = ['EC02','EC04','EE03','WC01','WW02']
new12normal = ['WC02','CC01']
new12noise = ['EC05','EE01','WW01','WW04','CC03','EC01','CC12']
new12other = ['CC10','EC03']
valid = ['CC01','CC09','CC10','CC12','EC05','EE01','WC02','WW01','WW04']
invalid = ['CC03','EC01','EC03']
old12 = old12normal + old12noise
mid6 = mid6noise + mid6other
new12 = new12normal + new12noise + new12other
laternoise = new12noise + mid6noise
driftsta1 = ['EC02','EC04','EE03','WC01','WW02']
driftsta2 = ['CC03','EC01','EC03']
#extrasta = ['EC01']
extrasta = []
    
def get_para(head='../paraP.cfg',ccch='Z'):
    # type: P, Z, PandZ, PorZ
    '''
    try:
        opts, args = getopt.getopt(sys.argv[1:], "C:P:")
        for op, value in opts:
            if op == '-C':
                ccch = value
            elif op == '-P':
                head = value
            else:
                Usage()
    except:
        pass
    # check channel        
    if channel != 'P' and channel != 'Z':
        print('No such channel: %s!'%channel)
        sys.exit(1)
    '''
    # check head file        
    if os.path.isfile(head):
        config.read(head)
    else:
        print("head file doesn't exist!")
        
    prefix = head.replace(os.path.basename(head),'')
    out_path = os.path.join(prefix,config.get('path', 'out_path'))
    image_path = os.path.join(out_path,config.get('path', 'image_dir'))
    path = {}
    path['P'] = os.path.join(out_path,config.get('path', 'P_dir'))
    path['Z'] = os.path.join(out_path,config.get('path', 'Z_dir'))
    logfile = os.path.join(prefix,config.get('path', 'logpath'),'evt_%s.log'%ccch)
    time_before = config.getfloat('para', 'time_beforetauP')
    time_after = config.getfloat('para', 'time_aftertauP')
    
    return ccch, path, image_path, logfile, time_before, time_after

def sortAB(A,B):
    B = [b for _,b in sorted(zip(A,B))]
    return B

def rms(data):
    return np.sqrt(np.mean(np.power(data,2)))

def SNR(data,fs,t_boundary,time_before,time_after):
    tb = t_boundary-time_before
    te = t_boundary+time_after
    nb = int(tb*fs-1)
    nc = int(t_boundary*fs-1)
    ne = int(te*fs-1)

    n_win = data[nb:nc]
    s_win = data[nc+1:ne+1]
    
    if rms(n_win) == 0: return 0

    return max(np.abs(s_win))/rms(n_win)

class Opts:
    def __init__(self):
        self.allevtlog = []
        self.evtlog = []
        self.ccwin = [-5, 5]
        self.evtnum = 0
        self.evtid = 0
        self.freq = [0.1, 2]
        self.noisegateP = 0
        self.noisegateZ = 0
        self.path = {}
        self.image_path = ""
        self.time_before = 50
        self.time_after = 50
        self.xlim = [-20, 20]
        self.dist_lim = [30, 95]
        self.cc_lagmax = 1.5
        self.amp_change_ratio = 1.4
        self.evtopts_lst = []
        self.ccch = ''
        self.mag_thre = 0
        self.norm = 0.4
        self.lines = []
        self.ccwinlst = []
        self.freqlst = []
        self.goodstach = []

class Evtopts:
    def __init__(self):
        self.staname = []
        self.channel = []
        self.dist = []
        self.trs = []
        self.trsfull = []
        self.trsbak = []
        self.norm = 1
        self.scale = 1
        self.trnum = 0
        self.npts = 0
        self.fs = 0
        self.snr = []
        self.goodtr = []
        self.waveform = []
        self.waveform_full = []
        self.dt = []
        self.dcstd = []
        self.dvcstd = []
        self.rsstd = []
        self.ylim = [0, 22]
        self.evla = 0
        self.evlo = 0
        self.evdp = 0
        self.mag = 0
        self.evttime = ''
        self.ccwin = []
        self.ccmax = []
        self.badtype = []
        self.freq = []

def submit(self,text,type):
    if type == 'cclag':
        self.opts.cc_lagmax = float(text)
    elif type == 'evtid':
        opts = self.opts
        evtid = int(text)
        if evtid < 0 or evtid >= opts.evtnum:
            print('No such event!')
            return
        opts.evtid = evtid
        self.evtread()
    elif type == 'ccwin':
        self.evtopts.ccwin = [float(tmp) for tmp in text.split(',')]
        self.cutwin()
        self.plotwave()
    elif type == 'freq':
        self.evtopts.freq = [float(tmp) for tmp in text.split(',')]
        self.changefreq()

#def submitcclag(self, text):
    #self.opts.cc_lagmax = float(text)

#def submitevtid(self,text):
    

class Plottrfig():
    
    fig = plt.figure(figsize=(12, 16), dpi=60)
    axreMccc = plt.axes([0.01, 0.91, 0.07, 0.03])
    axReloadEvt = plt.axes([0.22, 0.91, 0.08, 0.03])
    axsnrcull = plt.axes([0.32, 0.91, 0.06, 0.03])
    axampcull = plt.axes([0.40, 0.91, 0.06, 0.03])
    axdelete = plt.axes([0.48, 0.91, 0.1, 0.03])
    axfinish = plt.axes([0.6, 0.91, 0.07, 0.03])
    axback = plt.axes([0.7, 0.91, 0.07, 0.03])
    axnext = plt.axes([0.79, 0.91, 0.07, 0.03])
    #axnormalize = plt.axes([0.6, 0.91, 0.07, 0.03])
    #axcclag = plt.axes([0.17, 0.91, 0.03, 0.03])
    axfreq = plt.axes([0.14, 0.91, 0.07, 0.03])
    axccwin = plt.axes([0.04, 0.001, 0.07, 0.03])
    axevtid = plt.axes([0.9, 0.91, 0.07, 0.03])
    
    breMccc = Button(axreMccc, 'reMccc')
    bReloadEvt = Button(axReloadEvt, 'reloadEvt')
    bsnrcull = Button(axsnrcull, 'snrcull')
    bampcull = Button(axampcull, 'ampcull')
    bdelete = Button(axdelete, 'deleteBadTrs')
    bfinish = Button(axfinish, 'finish')
    bback = Button(axback, 'Back(n)')
    bnext = Button(axnext, 'Next(N)')
    
    #bnormalize = Button(axnormalize, 'normalize')
    
    '''
    axPCAcull = plt.axes([0.81, 0.91, 0.07, 0.03])
    axampcull = plt.axes([0.71, 0.91, 0.07, 0.03])
    axfinish = plt.axes([0.91, 0.91, 0.07, 0.03])
    axPlot = plt.axes([0.1, 0.91, 0.07, 0.03])
    axrePlot = plt.axes([0.2, 0.91, 0.07, 0.03])
    axdelete = plt.axes([0.3, 0.91, 0.14, 0.03])
    
    
    bPCAcull = Button(axPCAcull, 'PCAcull')
    bampcull = Button(axampcull, 'ampcull')
    bfinish = Button(axfinish, 'Finish')
    bplot = Button(axPlot, 'Plot RFs')
    breplot = Button(axrePlot, 'rePlot')
    bdelete = Button(axdelete, 'deletebadRFs')
    '''
    
    def __init__(self, opts):
        
        self.ax = plt.axes([0.2, 0.05, 0.15, 0.84])
        self.ax_full = plt.axes([0.4, 0.05, 0.3, 0.84])
        self.ax_dist = plt.axes([0.75, 0.05, 0.16, 0.84])
        self.opts = opts
        self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.fig.canvas.mpl_connect('key_press_event', self.keyboard)
        
        # buttons
        self.breMccc.on_clicked(self.mccc)
        self.bdelete.on_clicked(self.deletebadtrs)
        self.bReloadEvt.on_clicked(self.evtread)
        self.bback.on_clicked(self.Back)
        self.bnext.on_clicked(self.Next)
        self.bsnrcull.on_clicked(self.checksnr)
        self.bampcull.on_clicked(self.ampcull)
        self.bfinish.on_clicked(self.finish)
        #self.bnormalize.on_clicked(self.normalize)
        
        # text boxes
        self.txtfreq = TextBox(self.axfreq,'Freq[Hz]:',initial='%.1f,%.1f'%(opts.freq[0],opts.freq[1]))
        self.txtfreq.on_submit(lambda text: submit(self,text,type='freq'))
        #self.txtcclag = TextBox(self.axcclag,'maxlag[s]:',initial='%.1f'%opts.cc_lagmax)
        #self.txtcclag.on_submit(lambda text: submit(self,text,type='cclag'))
        self.txtccwin = TextBox(self.axccwin,'ccwin:',initial='%.1f,%.1f'%(opts.ccwin[0],opts.ccwin[1]))
        self.txtccwin.on_submit(lambda text: submit(self,text,type='ccwin'))
        self.txtevtid = TextBox(self.axevtid,'Evtid:',initial='')
        self.txtevtid.on_submit(lambda text: submit(self,text,type='evtid'))
        
        # plot waveform
        self.evtread()
        #self.ampcull()
        
    def onclick(self, event):
        #ax = self.ax
        ax = self.ax
        ax_full = self.ax_full
        evtopts = self.evtopts
        #print(event.ydata)
        
        if event.inaxes != ax and event.inaxes != ax_full:
            return
        click_idx = int(np.round(event.ydata))
        if click_idx > evtopts.trnum:
            return
        if evtopts.goodtr[click_idx-1] == 1:
            print("Selected "+evtopts.staname[click_idx-1])
            evtopts.goodtr[click_idx-1] = 0
            evtopts.waveform[click_idx-1].set_color('red')
            evtopts.waveform_full[click_idx-1].set_color('red')
            evtopts.badtype[click_idx-1] = 'click'
        else:
            print("Cancelled "+evtopts.staname[click_idx-1])
            evtopts.goodtr[click_idx-1] = 1
            evtopts.waveform[click_idx-1].set_color('black')
            evtopts.waveform_full[click_idx-1].set_color('black')
            evtopts.badtype[click_idx-1] = ''
        plt.draw()
        
    def keyboard(self, event):
        opts = self.opts
        evtopts = self.evtopts
        #print(event.key)
        if event.key == 'up':
            evtopts.scale *= opts.amp_change_ratio
            self.plotwave()
            
        elif event.key == 'down':
            evtopts.scale /= opts.amp_change_ratio
            self.plotwave()
            
        elif event.key == 'shift+up':
            evtopts.scale *= 5*opts.amp_change_ratio
            self.plotwave()
            
        elif event.key == 'shift+down':
            evtopts.scale /= 5*opts.amp_change_ratio
            self.plotwave()
            
        elif event.key == 'N':
            self.Lastpicked()
            
        elif event.key == 'n':
            self.Nextpicked()
            
        elif event.key == 's':
            self.save()
        
        #elif event.key == 'd':
        #    self.opts.lines[self.opts.evtid] = []
        #    self.opts.freqlst[self.opts.evtid] = np.zeros(2)
        #    self.opts.ccwinlst[self.opts.evtid] = np.zeros(2)
        #    self.opts.goodstach[self.opts.evtid] = []
        #    print('Delete evt%d'%self.opts.evtid)
        elif event.key == 'u':
            print('picked events:')
            for i in range(opts.evtnum):
                if opts.lines[i] != []: print('evt %d'%i)
        
        # frequency
        elif event.key == 'w':
            evtopts.freq = [0.1, 0.5]
            self.changefreq()
            
        elif event.key == 'e':
            evtopts.freq = [0.1, 1]
            self.changefreq()
        
        elif event.key == 'r':
            evtopts.freq = [0.1, 2]
            self.changefreq()
        
        elif event.key == 't':
            evtopts.freq = [0.5, 2]
            self.changefreq()
        
        elif event.key == 'y':
            evtopts.freq = [1, 2]
            self.changefreq()
        else:
            return
        
        
    def evtread(self, event=[]):
        opts = self.opts
        evtopts = Evtopts()
        evtopts.evttime = opts.evtlog[opts.evtid]
        if opts.ccch == 'Z' or opts.ccch == 'P':    
            saclst = os.path.join(opts.path[opts.ccch],evtopts.evttime,'%s*'%evtopts.evttime)
            st = obspy.read(saclst)
        elif opts.ccch == 'PorZ':
            st = []
            saclst_P = glob.glob(os.path.join(opts.path['P'],evtopts.evttime,'%s*'%evtopts.evttime))
            saclst_Z = glob.glob(os.path.join(opts.path['Z'],evtopts.evttime,'%s*'%evtopts.evttime))
            saclst = saclst_P + saclst_Z
            for sac in saclst: st.append(obspy.read(sac)[0])
            for tr in st:
                if tr.stats.channel == 'P': tr.data *= -1
        '''
        elif opts.ccch == 'PandZ':
            st = []
            evt_path_P = os.path.join(opts.path['P'],evtopts.evttime)
            evt_path_Z = os.path.join(opts.path['Z'],evtopts.evttime)
            sac_prefix_P = [os.path.basename(sac).remove('.P') for sac in glob.glob(os.path.join(evt_path_P,'%s*'%evtopts.evttime))]
            sac_prefix_Z = [os.path.basename(sac).remove('.Z') for sac in glob.glob(os.path.join(evt_path_Z,'%s*'%evtopts.evttime))]
            sac_prefix = sac_prefix_P[:]
            for prefix in sac_prefix_P:
                if not prefix in sac_prefix_Z:
                    sac_prefix.remove(prefix)
            for prefix in sac_prefix:
                trP = obspy.read(os.path.join(evt_path_P,prefix+'.P'))[0].copy()
                trZ = obspy.read(os.path.join(evt_path_Z,prefix+'.Z'))[0].copy()
                trZ.data = trZ.data - trP.data
                trZ.stats.channel = opts.ccch
                st.append(trZ)
        '''
        
        evtopts.trs = sorted(st, key=lambda x:x.stats.sac.gcarc)
        evtopts.trsbak = copy.deepcopy(evtopts.trs)
        evtopts.trnum = len(evtopts.trs)
        time = np.linspace(-opts.time_before,opts.time_after,num=evtopts.trs[0].stats.npts,endpoint=True)
        evtopts.fs = float(evtopts.trs[0].stats.sampling_rate)
        if opts.freqlst[opts.evtid].all():
            evtopts.freq = opts.freqlst[opts.evtid][:]
        else:
            evtopts.freq = opts.freq[:]
        if opts.ccwinlst[opts.evtid].all():
            evtopts.ccwin = opts.ccwinlst[opts.evtid][:]
        else:
            evtopts.ccwin = opts.ccwin[:]
        
        for i in range(evtopts.trnum):
            tr = evtopts.trs[i]
            tr.filter('bandpass',freqmin=evtopts.freq[0],freqmax=evtopts.freq[1],zerophase=True)
            evtopts.trsfull.append(tr.copy())
            evtopts.snr.append(SNR(tr.data,evtopts.fs,opts.time_before-5,50,15))
            tr.data = tr.data[np.where((time>=evtopts.ccwin[0])&(time<=evtopts.ccwin[1]))]
            tr.taper(max_percentage=0.05,type='hann')
            evtopts.staname.append(tr.stats.station)
            evtopts.channel.append(tr.stats.channel)
            evtopts.dist.append(tr.stats.sac.gcarc)
        evtopts.npts = evtopts.trs[0].stats.npts
        
        if opts.goodstach[opts.evtid] == []:
            evtopts.goodtr = np.ones(evtopts.trnum).astype('int')
        else:
            evtopts.goodtr = np.zeros(evtopts.trnum).astype('int')
            for i in range(evtopts.trnum):
                if evtopts.staname[i]+'.'+evtopts.channel[i] in opts.goodstach[opts.evtid]:
                    evtopts.goodtr[i] = 1
            
        evtopts.dt = np.zeros(evtopts.trnum)
        evtopts.dcstd = np.zeros(evtopts.trnum)
        evtopts.dvcstd = np.zeros(evtopts.trnum)
        evtopts.ccmax = np.zeros(evtopts.trnum)
        evtopts.rsstd = np.zeros(evtopts.trnum)
        evtopts.badtype = ['']*evtopts.trnum
        evtopts.ylim = [0, evtopts.trnum+3]
        evtopts.waveform = [[] for i in range(evtopts.trnum)]
        evtopts.waveform_full = [[] for i in range(evtopts.trnum)]
        evtopts.evla = evtopts.trs[0].stats.sac.evla
        evtopts.evlo = evtopts.trs[0].stats.sac.evlo
        evtopts.evdp = evtopts.trs[0].stats.sac.evdp
        evtopts.mag = evtopts.trs[0].stats.sac.mag
        self.evtopts = evtopts
        self.normalize()
        self.plotwave()
        #self.plotwave()
        #self.checksnr()
        #self.ampcull()
        #self.mccc()
    
    def trspectrogram(self):
        opts = self.opts
        evtopts = self.evtopts
        fig_spec = plt.figure(figsize=(16,16))
        for i in range(evtopts.trnum):
            ax_spec = fig_spec.add_subplot(np.ceil(evtopts.trnum/3.),3,i+1)
            f,t,psd = spectrogram(evtopts.trsbak[i].data,int(evtopts.fs),nperseg=int(10*evtopts.fs),noverlap=int(9*evtopts.fs))
            psd = 10*np.log10(psd)
            ax_spec.semilogy()
            levels = np.linspace(0.4*psd.min(),psd.max(),50)
            cs = ax_spec.contourf(t-opts.time_before,f,psd,levels,cmap=plt.cm.rainbow)
            cbar = fig_spec.colorbar(cs,ax=ax_spec)
            cbar.ax.set_ylabel('PSD[dB]')
            ax_spec.set(xlim=opts.xlim,ylim=[0.1,2],xlabel='t[s]',ylabel='freq[Hz]',title=evtopts.staname[i]+'.'+evtopts.channel[i])
        plt.draw()
        #plt.show()
        #plt.close(fig_spec)
    
    def changefreq(self):
        opts = self.opts
        evtopts = self.evtopts
        time = np.linspace(-opts.time_before,opts.time_after,num=evtopts.trsfull[0].stats.npts,endpoint=True)
        evtopts.trs = []
        evtopts.trsfull = []
        evtopts.snr = []
        #evtopts.goodtr = np.ones(evtopts.trnum).astype('int')
        for i in range(evtopts.trnum):
            tr = evtopts.trsbak[i].copy()
            tr.filter('bandpass',freqmin=evtopts.freq[0],freqmax=evtopts.freq[1],zerophase=True)
            evtopts.snr.append(SNR(tr.data,evtopts.fs,opts.time_before-5,50,15))
            evtopts.trsfull.append(tr.copy())
            tr.data = tr.data[np.where((time>=evtopts.ccwin[0])&(time<=evtopts.ccwin[1]))]
            tr.taper(max_percentage=0.05,type='hann')
            evtopts.trs.append(tr)
        print('Adjust frequency: %.2f-%.2f[Hz]'%(evtopts.freq[0],evtopts.freq[1]))
        self.normalize()
        self.plotwave()
    
    def cutwin(self):
        opts = self.opts
        evtopts = self.evtopts
        time = np.linspace(-opts.time_before,opts.time_after,num=evtopts.trsfull[0].stats.npts,endpoint=True)
        for i in range(evtopts.trnum):
            tr = evtopts.trs[i]
            tr.data = evtopts.trsfull[i].data.copy()
            tr.data = tr.data[np.where((time>=evtopts.ccwin[0])&(time<=evtopts.ccwin[1]))]
            tr.taper(max_percentage=0.05,type='hann')
        evtopts.npts = evtopts.trs[0].stats.npts
    
    def normalize(self, event=[]):
        opts = self.opts
        evtopts = self.evtopts
        good_id = np.where(evtopts.goodtr==1)[0]
        if len(good_id) == 0: good_id = np.arange(0,evtopts.trnum)
        if len(np.unique(evtopts.channel)) == 1:
            absampmax = []
            for i in good_id:
                absampmax.append(np.abs(evtopts.trs[i].data).max())
            evtopts.norm = 1./np.median(absampmax)*opts.norm
            
        elif len(np.unique(evtopts.channel)) == 2:
            channel = [evtopts.channel[i] for i in good_id]
            trs = [evtopts.trs[i] for i in good_id]
            if len(np.unique(channel)) == 2:
                ch1, ch2 = np.unique(channel)
                if evtopts.channel.count(ch1) < evtopts.channel.count(ch2): ch1, ch2 = ch2, ch1
                absampmax1, absampmax2 = [], []
                for i in range(len(channel)):
                    if channel[i] == ch1: absampmax1.append(np.abs(trs[i].data).max())
                    elif channel[i] == ch2: absampmax2.append(np.abs(trs[i].data).max())
                evtopts.norm = 1./np.median(absampmax1)*opts.norm
                norm_2to1 = np.median(absampmax1)/np.median(absampmax2)
                for i in range(evtopts.trnum):
                    if evtopts.channel[i] == ch2:
                        evtopts.trs[i].data *= norm_2to1
                        evtopts.trsfull[i].data *= norm_2to1
                    
            elif len(np.unique(channel)) == 1:
                ch1 = np.unique(channel)[0]
                absampmax1, absampmax2 = [], []
                for i in range(len(channel)): absampmax1.append(np.abs(trs[i].data).max())
                evtopts.norm = 1./np.median(absampmax1)*opts.norm
                for i in range(evtopts.trnum):
                    if evtopts.channel[i] != ch1: absampmax2.append(np.abs(evtopts.trs[i].data).max())
                norm_2to1 = np.median(absampmax1)/np.median(absampmax2)
                for i in range(evtopts.trnum):
                    if evtopts.channel[i] != ch1:
                        evtopts.trs[i].data *= norm_2to1
                        evtopts.trsfull[i].data *= norm_2to1
                
    def mccc(self, event=[]):
        #self.normalize()
        print('cc_lagmax = %.1f[s]'%self.opts.cc_lagmax)
        opts = self.opts
        evtopts = self.evtopts
        good_id = np.where(evtopts.goodtr==1)[0]
        bad_id = np.where(evtopts.goodtr==0)[0]
        matrice = np.zeros([evtopts.npts,len(good_id)])
        if len(good_id) < 3:
            print("less than 3 traces, can't do mccc")
            return
        for k in range(len(good_id)):
            matrice[:,k] = evtopts.trs[good_id[k]].data
        #matrice = matlab.double(matrice.tolist())
        dcor,dcstd,dvcstd,_,rsstd = cc.mccc(matrice,1./evtopts.fs,lagmax=opts.cc_lagmax,ifwt=True)
        
        #dcor = np.array(dcor).reshape(1,-1)[0] # smaller - phase earlier
        #dcstd = np.array(dcstd).reshape(1,-1)[0]
        #dvcstd = np.array(dvcstd).reshape(1,-1)[0]
        print('maxdiff = %.4f'%(dcor.max()-dcor.min()))
        print('dcor:')
        for k in range(len(good_id)):
            #print('%s %.3f'%(evtopts.staname[good_id[k]],dcor[k]))
            evtopts.dt[good_id[k]] = dcor[k]
            evtopts.dcstd[good_id[k]] = dcstd[k]
            evtopts.dvcstd[good_id[k]] = dvcstd[k]
            evtopts.rsstd[good_id[k]] = rsstd[k]
        for k in bad_id:
            evtopts.dt[k] = 0
            evtopts.dcstd[k] = 0
            evtopts.dvcstd[k] = 0
            evtopts.rsstd[k] = 0
        self.ccmax_cal()
        self.plotwave()
    
    def ccmax_cal(self):
        opts = self.opts
        evtopts = self.evtopts
        stack = np.zeros(evtopts.npts)
        time_axis = np.linspace(evtopts.ccwin[0], evtopts.ccwin[1], evtopts.npts)
        good_id = np.where(evtopts.goodtr==1)[0]
        for i in good_id:
            stack += np.interp(time_axis,time_axis-evtopts.dt[i],evtopts.trs[i].data)
        stack /= len(good_id)
        for i in range(evtopts.trnum):
            #cc = correlate(evtopts.trs[i].data,stack,shift=evtopts.npts-1,demean=False,normalize='naive')
            #shift, evtopts.ccmax[i] = xcorr_max(cc,abs_max=False)
            c, shift = cc.xcorr(evtopts.trs[i].data,stack,lagmax=int(evtopts.fs*opts.cc_lagmax),scale='coeff')
            pk = np.where(c==max(c))[0]
            evtopts.ccmax[i] = max(c)
            shift = shift[pk]
            if evtopts.goodtr[i] == 1:
                print('%s: dt=%.3f[s], std=%.3f, ccmax=%.3f, shift=%.3f'%(evtopts.staname[i],evtopts.dt[i],evtopts.dvcstd[i],evtopts.ccmax[i],shift*1./evtopts.fs))
    
    def checksnr(self, event=[]):
        opts = self.opts
        evtopts = self.evtopts
        for i in range(evtopts.trnum):
            if evtopts.goodtr[i] == 1: 
                if evtopts.channel[i] == 'P' and evtopts.snr[i] < opts.noisegateP:
                    evtopts.goodtr[i] = 0
                    evtopts.waveform[i].set_color('red')
                    evtopts.waveform_full[i].set_color('red')
                    print('bad trace of station %s.%s : snr %.2f is smaller than noisegate %.2f'\
                          %(evtopts.staname[i],evtopts.channel[i],evtopts.snr[i],opts.noisegateP))
                    evtopts.badtype[i] = 'snr'
                if evtopts.channel[i] == 'Z' and evtopts.snr[i] < opts.noisegateZ:
                    evtopts.goodtr[i] = 0
                    evtopts.waveform[i].set_color('red')
                    evtopts.waveform_full[i].set_color('red')
                    print('bad trace of station %s.%s : snr %.2f is smaller than noisegate %.2f'\
                          %(evtopts.staname[i],evtopts.channel[i],evtopts.snr[i],opts.noisegateZ))
                    evtopts.badtype[i] = 'snr'
        plt.draw()
        
    def cull(self, event=[]):
        self.checksnr()
        self.ampcull()
    
    def ampcull(self, event=[]):
        self.normalize()
        opts = self.opts
        evtopts = self.evtopts
        for i in range(evtopts.trnum):
            if evtopts.goodtr[i] == 1 or evtopts.goodtr[i] == 0:
                tr = evtopts.trs[i]
                maxabsamp = np.abs(tr.data).max()
                if maxabsamp == 0:
                    evtopts.goodtr[i] = 0
                    evtopts.waveform[i].set_color('red')
                    evtopts.waveform_full[i].set_color('red')
                    print('bad trace of station %s.%s : constant waveform'%(evtopts.staname[i],evtopts.channel[i]))
                    evtopts.badtype[i] = 'constant'
                elif maxabsamp < 0.1*1./evtopts.norm*opts.norm:
                    evtopts.goodtr[i] = 0
                    evtopts.waveform[i].set_color('red')
                    evtopts.waveform_full[i].set_color('red')
                    print('bad trace of station %s.%s : too small amp'%(evtopts.staname[i],evtopts.channel[i]))
                    evtopts.badtype[i] = 'loamp'
                elif maxabsamp > 10*1./evtopts.norm*opts.norm:
                    evtopts.goodtr[i] = 0
                    evtopts.waveform[i].set_color('red')
                    evtopts.waveform_full[i].set_color('red')
                    print('bad trace of station %s.%s : too large amp'%(evtopts.staname[i],evtopts.channel[i]))
                    evtopts.badtype[i] = 'hiamp'
        #self.normalize()
        #print([(t1,t2,t3) for t1,t2,t3 in zip(evtopts.goodtr,evtopts.staname,evtopts.snr)])
        #self.deletebadtrs(type='hiamp')
        #print([(t1,t2,t3) for t1,t2,t3 in zip(self.evtopts.goodtr,self.evtopts.staname,self.evtopts.snr)])
    
    def Next(self, event=[]):
        opts = self.opts
        #evtopts = self.evtopts
        if opts.evtid+1 >= opts.evtnum:
            print('Already the last event')
            return
        #evtopts.waveform = [[]]*evtopts.trnum
        #opts.evtopts_lst[opts.evtid] = copy.deepcopy(evtopts)
        opts.evtid += 1
        self.evtread()
        #self.ampcull()
        print('next event')    
        
    def Back(self, event=[]):
        opts = self.opts
        if opts.evtid-1 < 0:
            print('Already the first event')
            return
        opts.evtid -= 1
        self.evtread()
        #self.ampcull()
        #self.evtopts = opts.evtopts_lst[opts.evtid]
        #self.plotwave()
        print('last event')
        
    def Nextpicked(self,event=[]):
        opts = self.opts
        opts.evtid = (opts.evtid+1)%opts.evtnum
        while(opts.lines[opts.evtid]==[]):
            opts.evtid = (opts.evtid+1)%opts.evtnum
        self.evtread()
        print('next picked event')
        
    def Lastpicked(self,event=[]):
        opts = self.opts
        opts.evtid = (opts.evtid-1)%opts.evtnum
        while(opts.lines[opts.evtid]==[]):
            opts.evtid = (opts.evtid-1)%opts.evtnum
        self.evtread()
        print('last picked event')
                
    def plotdist(self):
        ax_dist = self.ax_dist
        opts = self.opts
        evtopts = self.evtopts
        ax_dist.cla()
        ax_dist.set_xlabel("Distance[deg]")
        ax_dist.scatter(evtopts.dist, np.arange(evtopts.trnum)+1)
        ax_dist.set_xlim(opts.dist_lim)
        ax_dist.set_xticks(np.arange(opts.dist_lim[0],opts.dist_lim[1],10))
        ax_dist.set_ylim(evtopts.ylim)
        ax_dist.set_yticks(np.arange(evtopts.trnum)+1)
        ax_dist.set_yticklabels([])
        ax_dist.grid()
    
    def plotwave(self, event=[]):
        ax = self.ax
        ax_full = self.ax_full
        opts = self.opts
        evtopts = self.evtopts
        
        ax.cla()
        ax.set_xlim(np.array(evtopts.ccwin)+np.array([-3,3]))
        ax.set_ylim(evtopts.ylim)
        ax.set_ylabel("sta", fontsize=20)
        ax.set_xlabel("Time[s]")
        ax.set_title("cc segment")
        ax.grid()
        
        ax_full.cla()
        ax_full.set_xlim(opts.xlimfull)
        ax_full.set_ylim(evtopts.ylim)
        #ax_full.set_ylabel("sta", fontsize=20)
        ax_full.set_xlabel("Time[s]")
        ax_full.set_title("Traces (%.2f-%.2fHz)"%(evtopts.freq[0],evtopts.freq[1]))
        ax_full.grid()
        
        trs = [[], []]
        trsfull = [[], []]
        trsbak = [[], []]
        goodtr = [[], []]
        staname = [[], []]
        channel = [[], []]
        snr = [[], []]
        dist = [[], []]
        dt = [[], []]
        dcstd = [[], []]
        dvcstd = [[], []]       
        ccmax = [[], []]
        rsstd = [[], []]
        badtype = [[], []]
        sortvalue = [[], []]
        
        # badtf - front, goodtf - back
        for i in range(evtopts.trnum):
            trs[evtopts.goodtr[i]].append(evtopts.trs[i])
            trsfull[evtopts.goodtr[i]].append(evtopts.trsfull[i])
            trsbak[evtopts.goodtr[i]].append(evtopts.trsbak[i])
            staname[evtopts.goodtr[i]].append(evtopts.staname[i])
            channel[evtopts.goodtr[i]].append(evtopts.channel[i])
            snr[evtopts.goodtr[i]].append(evtopts.snr[i])
            dist[evtopts.goodtr[i]].append(evtopts.dist[i])
            dt[evtopts.goodtr[i]].append(evtopts.dt[i])
            dcstd[evtopts.goodtr[i]].append(evtopts.dcstd[i])
            dvcstd[evtopts.goodtr[i]].append(evtopts.dvcstd[i])
            ccmax[evtopts.goodtr[i]].append(evtopts.ccmax[i])
            rsstd[evtopts.goodtr[i]].append(evtopts.rsstd[i])
            badtype[evtopts.goodtr[i]].append(evtopts.badtype[i])
            goodtr[evtopts.goodtr[i]].append(evtopts.goodtr[i])
            #if len(np.unique(evtopts.ccmax)) > 1:
            sortvalue[evtopts.goodtr[i]].append('%.5f'%evtopts.ccmax[i]+evtopts.channel[i]+evtopts.staname[i])
            #else:
            #    sortvalue[evtopts.goodtr[i]].append(evtopts.channel[i]+evtopts.staname[i][:-1])
            
        for tmp in [trs,trsfull,trsbak,staname,channel,dist,snr,goodtr,dt,dcstd,dvcstd,ccmax,rsstd,badtype]:
            tmp[0] = sortAB(sortvalue[0][:],tmp[0])
            tmp[1] = sortAB(sortvalue[1][:],tmp[1])
        evtopts.trs = trs[0] + trs[1]
        evtopts.trsfull = trsfull[0] + trsfull[1]
        evtopts.trsbak = trsbak[0] + trsbak[1]
        evtopts.staname = staname[0] + staname[1]
        evtopts.channel = channel[0] + channel[1]
        evtopts.snr = snr[0] + snr[1]
        evtopts.goodtr = goodtr[0] + goodtr[1]
        evtopts.goodtr = np.array(evtopts.goodtr).astype('int')
        evtopts.dist = dist[0] + dist[1]
        evtopts.dt = dt[0] + dt[1]
        evtopts.dcstd = dcstd[0] + dcstd[1]
        evtopts.dvcstd = dvcstd[0] + dvcstd[1]
        evtopts.ccmax = ccmax[0] + ccmax[1]
        evtopts.rsstd = rsstd[0] + rsstd[1]
        evtopts.badtype = badtype[0] + badtype[1]
        
        time_axis = np.linspace(evtopts.ccwin[0], evtopts.ccwin[1], evtopts.npts)
        time_axis_full = np.linspace(-opts.time_before,opts.time_after,num=evtopts.trsfull[0].stats.npts,endpoint=True)
        stack = np.zeros(evtopts.npts)
        stack_full = np.zeros(evtopts.trsfull[0].stats.npts)
        
        for i in range(evtopts.trnum):
            tr = evtopts.trs[i]
            trfull = evtopts.trsfull[i]
            #amp_axis = rf.data*opts.enf+i+1
            #amp_axis = rf.data/1.2/np.max(rf.data)+i+1
            amp_axis = tr.data*evtopts.norm*evtopts.scale+i+1
            amp_axis_full = trfull.data*evtopts.norm*evtopts.scale+i+1
            
            if evtopts.goodtr[i] == 0:
                evtopts.waveform[i] = ax.plot(time_axis-evtopts.dt[i], amp_axis, color="red", linewidth=1)[0]
                evtopts.waveform_full[i] = ax_full.plot(time_axis_full, amp_axis_full, color="red", linewidth=1)[0]
            elif evtopts.goodtr[i] == 1:
                if evtopts.staname[i] not in extrasta:
                    evtopts.waveform[i] = ax.plot(time_axis-evtopts.dt[i], amp_axis, color="black", linewidth=1)[0]
                    evtopts.waveform_full[i] = ax_full.plot(time_axis_full, amp_axis_full, color="black", linewidth=1)[0]
                else:
                    evtopts.waveform[i] = ax.plot(time_axis-evtopts.dt[i], amp_axis, color="blue", linewidth=1)[0]
                    evtopts.waveform_full[i] = ax_full.plot(time_axis_full, amp_axis_full, color="blue", linewidth=1)[0]
                stack += np.interp(time_axis,time_axis-evtopts.dt[i],tr.data)
                stack_full += np.interp(time_axis_full,time_axis_full-evtopts.dt[i],trfull.data)
        
        
        stack /= len(np.where(evtopts.goodtr==1)[0])
        stack_full /= len(np.where(evtopts.goodtr==1)[0])
        ax.plot(time_axis,stack*evtopts.norm*evtopts.scale+evtopts.trnum+2,color="black",linewidth=2)
        
        ax_full.plot(time_axis_full,stack_full*evtopts.norm*evtopts.scale+evtopts.trnum+2,color="black",linewidth=2)
        ax_full.vlines(evtopts.ccwin[0],evtopts.ylim[0],evtopts.ylim[1]+3,linewidth=1.5,linestyle='dashed',color='blue')
        ax_full.vlines(evtopts.ccwin[1],evtopts.ylim[0],evtopts.ylim[1]+3,linewidth=1.5,linestyle='dashed',color='blue')
        
        add_ytick = []
        add_ylabel = []
        try:
            ref_st = obspy.read('TAOE/%s_TAOE.BHZ'%evtopts.evttime)[0].copy()
            add_ytick.append(evtopts.trnum+3)
            add_ylabel.append('TAOE') 
            ref_st.filter('bandpass',freqmin=evtopts.freq[0],freqmax=evtopts.freq[1],zerophase=True)
            ref_st.taper(max_percentage=0.05,type='hann')
            time_axis_full = np.linspace(-opts.time_before,opts.time_after,num=ref_st.stats.npts,endpoint=True)
            if opts.ccch != 'P':
                ax_full.plot(time_axis_full,ref_st.data/np.abs(ref_st.data).max()*opts.norm*evtopts.scale+evtopts.trnum+3,color="black",linewidth=2)
            else:
                ax_full.plot(time_axis_full,-ref_st.data/np.abs(ref_st.data).max()*opts.norm*evtopts.scale+evtopts.trnum+3,color="black",linewidth=2)
            #ax_full.plot(time_axis_full,ref_st.data/np.abs(ref_st.data).max()*opts.norm+evtopts.trnum+3,color="black",linewidth=2)
            ref_st.data = ref_st.data[np.where((time_axis_full>=evtopts.ccwin[0])&(time_axis_full<=evtopts.ccwin[1]))]
            ref_st.taper(max_percentage=0.05,type='hann')
            time_axis = np.linspace(evtopts.ccwin[0], evtopts.ccwin[1], num=ref_st.stats.npts,endpoint=True)
            if opts.ccch != 'P':
                ax.plot(time_axis,ref_st.data/np.abs(ref_st.data).max()*opts.norm*evtopts.scale+evtopts.trnum+3,color="black",linewidth=2)
            else:
                ax.plot(time_axis,-ref_st.data/np.abs(ref_st.data).max()*opts.norm*evtopts.scale+evtopts.trnum+3,color="black",linewidth=2)
        except:
            print('No land station record')
        
        try:
            ref_st = obspy.read('synthetic/%s_syn.BHZ'%evtopts.evttime)[0].copy()
            if evtopts.trnum+3 in add_ytick:
                add_ytick.append(evtopts.trnum+4)
            else:
                add_ytick.append(evtopts.trnum+3)
            add_ylabel.append('synthetic') 
            ref_st.filter('bandpass',freqmin=evtopts.freq[0],freqmax=evtopts.freq[1],zerophase=True)
            ref_st.taper(max_percentage=0.05,type='hann')
            time_axis_full = np.linspace(-opts.time_before,opts.time_after,num=ref_st.stats.npts,endpoint=True)
            if opts.ccch != 'P':
                ax_full.plot(time_axis_full,ref_st.data/np.abs(ref_st.data).max()*opts.norm*evtopts.scale+add_ytick[-1],color="black",linewidth=2)
            else:
                ax_full.plot(time_axis_full,-ref_st.data/np.abs(ref_st.data).max()*opts.norm*evtopts.scale+add_ytick[-1],color="black",linewidth=2)
            #ax_full.plot(time_axis_full,ref_st.data/np.abs(ref_st.data).max()*opts.norm+evtopts.trnum+3,color="black",linewidth=2)
            ref_st.data = ref_st.data[np.where((time_axis_full>=evtopts.ccwin[0])&(time_axis_full<=evtopts.ccwin[1]))]
            ref_st.taper(max_percentage=0.05,type='hann')
            time_axis = np.linspace(evtopts.ccwin[0], evtopts.ccwin[1], num=ref_st.stats.npts,endpoint=True)
            if opts.ccch != 'P':
                ax.plot(time_axis,ref_st.data/np.abs(ref_st.data).max()*opts.norm*evtopts.scale+add_ytick[-1],color="black",linewidth=2)
            else:
                ax.plot(time_axis,-ref_st.data/np.abs(ref_st.data).max()*opts.norm*evtopts.scale+add_ytick[-1],color="black",linewidth=2)
        except:
            print('No synthetic')
            
        ax.set_yticks(np.arange(1,evtopts.trnum+1).tolist()+[evtopts.trnum+2]+add_ytick)
        ax.set_yticklabels([tmp1+'.'+tmp2 for tmp1,tmp2 in zip(evtopts.staname,evtopts.channel)]+['stack']+add_ylabel)
        ax_full.set_yticks(np.arange(1,evtopts.trnum+1).tolist()+[evtopts.trnum+2]+add_ytick)
        ylimax = int(np.max([evtopts.trnum+2]+add_ytick))
        ax.set_ylim([0,ylimax+1])
        ax_full.set_ylim([0,ylimax+1])
            
        if len(np.unique(evtopts.ccmax)) > 1:
            ax_full.set_yticklabels(['%.3f'%tmp for tmp in evtopts.ccmax])
        else:
            ax_full.set_yticklabels([])
        
        self.plotdist()
        #print(len(opts.lines[opts.evtid]))
        if opts.lines[opts.evtid] == []:
            ispicked = ''
        else:
            ispicked = '(Picked) '
        self.fig.suptitle("%sEvt%d/%d %s (Evdp: %5.2fkm, Mw: %.1f)" % (ispicked, opts.evtid, opts.evtnum-1, evtopts.evttime, evtopts.evdp, evtopts.mag), fontsize=20)
        #self.trspectrogram()

    def deletebadtrs(self, event=[], type='all'):
        # type : all, amp, snr, click
        self.evtopts.waveform = [[] for i in range(self.evtopts.trnum)]
        self.evtopts.waveform_full = [[] for i in range(self.evtopts.trnum)]
        evtopts = copy.deepcopy(self.evtopts)
        if type == 'all':
            bad_id = np.where(evtopts.goodtr==0)[0]
        else:
            badtype = np.array(evtopts.badtype)
            bad_id = np.where(badtype==type)[0]
        evtopts.trs = []
        evtopts.trsfull = []
        evtopts.trsbak = []
        evtopts.staname = []
        evtopts.channel = []
        evtopts.snr = []
        evtopts.dt = []
        evtopts.dcstd = []
        evtopts.dvcstd = []
        evtopts.dist = []
        evtopts.goodtr = []
        evtopts.ccmax = []
        evtopts.badtype = []
        
        for i in range(evtopts.trnum):
            if i in bad_id:
                print("Reject station "+self.evtopts.staname[i]+'.'+self.evtopts.channel[i])
                continue
            evtopts.trs.append(self.evtopts.trs[i])
            evtopts.trsfull.append(self.evtopts.trsfull[i])
            evtopts.trsbak.append(self.evtopts.trsbak[i])
            evtopts.staname.append(self.evtopts.staname[i])
            evtopts.channel.append(self.evtopts.channel[i])
            evtopts.snr.append(self.evtopts.snr[i])
            evtopts.dt.append(self.evtopts.dt[i])
            evtopts.dcstd.append(self.evtopts.dcstd[i])
            evtopts.dvcstd.append(self.evtopts.dvcstd[i])
            evtopts.dist.append(self.evtopts.dist[i])
            evtopts.goodtr.append(self.evtopts.goodtr[i])
            evtopts.ccmax.append(self.evtopts.ccmax[i])
            evtopts.badtype.append(self.evtopts.badtype[i])
        evtopts.goodtr = np.array(evtopts.goodtr).astype('int')
        evtopts.trnum = len(evtopts.trs)
        evtopts.waveform = [[] for i in range(evtopts.trnum)]
        evtopts.waveform_full = [[] for i in range(self.evtopts.trnum)]
        evtopts.ylim = [0, evtopts.trnum+3]
        self.evtopts = evtopts
        self.plotwave()
    
    def save(self, event=[]):
        opts = self.opts
        evtopts = self.evtopts
        good_id = np.where(evtopts.goodtr==1)[0]
        opts.lines[opts.evtid] = []
        opts.ccwinlst[opts.evtid] = np.zeros(2)
        opts.freqlst[opts.evtid] = np.zeros(2)
        opts.goodstach[opts.evtid] = []
        cf =  np.sqrt(evtopts.freq[0]*evtopts.freq[1])
        orid = opts.allevtlog.index(evtopts.evttime)
        for i in good_id:
            tr = evtopts.trs[i]
            rayp_rad = np.deg2rad(tr.stats.sac.user4*6371)
            gcarc = tr.stats.sac.gcarc
            baz = tr.stats.sac.baz
            sta = evtopts.staname[i]
            ch = evtopts.channel[i]
            nwk = tr.stats.network
            evla = tr.stats.sac.evla
            evlo = tr.stats.sac.evlo
            evdp = tr.stats.sac.evdp
            dt = evtopts.dt[i]
            dcstd = evtopts.dcstd[i]
            dvcstd = evtopts.dvcstd[i]
            acor = evtopts.ccmax[i]
            rsstd = evtopts.rsstd[i]
            opts.lines[opts.evtid].append('%16s %8.4f %8.3f %8.2f %8s %8s %8d %8.4f %8.4f %8.2f %8.3f %8.4f %8.4f %8.4f %8.5f %8.4f %8s\n'\
                      %(evtopts.evttime,rayp_rad,gcarc,baz,sta,nwk,orid,evla,evlo,evdp,dt,dcstd,dvcstd,acor,rsstd,cf,ch))
            opts.goodstach[opts.evtid].append(sta+'.'+ch)
        print('evt%d saved'%opts.evtid)
        print('cf:%.2f'%cf)
        #plt.savefig(os.path.join(opts.dirname,'%s_%s.png'%(evtopts.evttime,opts.ccch)))
        opts.ccwinlst[opts.evtid] = np.array(evtopts.ccwin)
        opts.freqlst[opts.evtid] = np.array(evtopts.freq)
        
    def finish(self, event=[]):
        opts = self.opts
        lines = []
        for i in range(opts.evtnum):
            lines += opts.lines[i]
        if os.path.exists(opts.filename):
            os.rename(opts.filename,os.path.join(opts.dirname,'%s_mag_%.1f_%.1f.bak'%(opts.ccch,opts.mag_thre[0],opts.mag_thre[1])))
        with open(opts.filename,'w+') as f:
            f.writelines(lines)
            f.close()
        print('Results written into %s'%opts.filename)
        np.savetxt(opts.winfreqfile,np.hstack((opts.ccwinlst,opts.freqlst)))

def main():
    opts = Opts()
    opts.ccch, opts.path, opts.image_path, logfile, opts.time_before, opts.time_after = get_para()
    opts.noisegateZ = 3
    opts.noisegateP = 3
    opts.cc_lagmax = 1
    opts.mag_thre = [6, 10]
    #opts.path = path['Z']
    evtlog = np.loadtxt(logfile,usecols=0,dtype='str')
    maglog = np.loadtxt(logfile,usecols=4)
    opts.allevtlog = evtlog.tolist()
    opts.evtlog = evtlog[np.where((maglog>=opts.mag_thre[0])&(maglog<opts.mag_thre[1]))].tolist()
    opts.evtnum = len(opts.evtlog)
    opts.evtid = 0 #PorZ 16
    opts.ccwin = [-3, 5]
    opts.xlim = [-8, 5]
    opts.xlimfull = [-20, 20]
    opts.dist_lim = [30, 95]
    opts.freq = [0.3, 0.6]
    opts.amp_change_ratio = 1.4
    opts.evtopts_lst = [[] for i in range(opts.evtnum)]
    opts.lines = [[] for i in range(opts.evtnum)]
    opts.ccwinlst = np.zeros([opts.evtnum,2])
    opts.freqlst = np.zeros([opts.evtnum,2])
    opts.dirname = 'cclog/%s_mag_%.1f_%.1f'%(opts.ccch,opts.mag_thre[0],opts.mag_thre[1])
    opts.filename = os.path.join(opts.dirname,'%s_mag_%.1f_%.1f.dat'%(opts.ccch,opts.mag_thre[0],opts.mag_thre[1]))
    opts.winfreqfile = os.path.join(opts.dirname,'%s_mag_%.1f_%.1f_winfreq.dat'%(opts.ccch,opts.mag_thre[0],opts.mag_thre[1]))
    opts.goodstach = [[] for i in range(opts.evtnum)]
    
    if os.path.exists(opts.filename):
        with open(opts.filename,'r') as f:
            for line in f.readlines():
                evtid = int(line.split()[6])
                opts.lines[opts.evtlog.index(opts.allevtlog[evtid])].append(line)
                line = line.replace('\n','').split()
                opts.goodstach[opts.evtlog.index(opts.allevtlog[evtid])].append(line[4]+'.'+line[-1])
                opts.goodstach[opts.evtlog.index(opts.allevtlog[evtid])] += [tmp+'.'+opts.ccch for tmp in extrasta]
            f.close()
        print('load %s'%opts.filename)
    
    if os.path.exists(opts.winfreqfile):
        opts.ccwinlst = np.loadtxt(opts.winfreqfile,usecols=(0,1))
        opts.freqlst = np.loadtxt(opts.winfreqfile,usecols=(2,3))
        print('load %s'%opts.winfreqfile)
    if not os.path.exists(opts.dirname): os.makedirs(opts.dirname)
    Plottrfig(opts)

if __name__ == "__main__":
    main()
    plt.show()
