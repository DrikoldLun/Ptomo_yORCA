ofile = '/home/lun/Desktop/Ptomo/tomo/BWTOMOG_atten_Vp_Vs_joint/Ltest_dT_PZmixwt.mat';
datfile = 'yORCA/data/PZmix.dat'
load(ofile)
Ltest.vr = Ltest.wvr;
resid2rough = 0.01:0.01:0.1;
semb = zeros(length(resid2rough),1);
for i = 1:length(resid2rough)
penalty = (Ltest.norm) + resid2rough(i)*(100-Ltest.vr);
[minA,x,y] = mingrid(penalty);
semb(i) = synth_test_ORCA(datfile,Ltest.damp(y),Ltest.smooth(x));
end