clear par 
% project details
global dbdir
dbdir = '/media/lun/easystore/UCSB/research/PacificORCA/Ptomo_new/tomo/BWTOMOG_atten_Vp_Vs_joint/yORCA/'; % include final slash
addpath('../matguts')
addpath('../')

if ~exist('par','var')==1, par=struct([]); end

par(1).phasecomp = 'PZ'; % 'PZ' or 'ST' or 'Sall' or 'BOTH_PZ_Sall'
% par(1).phasecomp = 'Sall'; % 'PZ' or 'ST' or 'Sall' or 'BOTH_PZ_Sall'
par.dattype = 'dT'; % 'dT' or 'dtstar'

%% results directories
% run([dbdir,dbname,'_startup.m']);
par.tomodatdir = [dbdir,'data/'];
par.figdir = [dbdir,'figs/'];
par.tomoresdir = [dbdir,'results/'];

[datfile,stafile,par.phases] = parse_phasecomp(par.tomodatdir,par.phasecomp,par.dattype);
datfile = 'yORCA/data/PZmix.dat'; %2,5
%datfile = 'yORCA/data/Z_mag_5.5_10.0.dat'; %1,5
%datfile = 'yORCA/data/PZ_direct_num.dat';%2,5
%datfile = 'yORCA/data/P_mag_6.0_10.0.dat'; %2,5
%datfile = 'yORCA/data/PZ_direct_both.dat'; %2,5

crustfile  = sprintf('%sstations_crust.dat',par.tomodatdir);

%% ---------------------------------
par.PS = find(strcmp({'P','S','B'},par.phasecomp(1))); % 1 for P, 2 for S, 3 for both
par.t_ts = find(strcmp({'dT','dtstar'},par.dattype)); % 1 for dT - v_inversion, 2 for dtstar - q_inversion

%% inversion running parms - zero for speed
REBUILD             = 0;
par.redo_crust_corr = REBUILD; % <== REDO IF DATA/MODEL CHANGES
par.build_K         = REBUILD; % <== REDO IF DATA/MODEL CHANGES
par.build_smooth    = 1; % <== REDO IF DATA/MODEL CHANGES
par.solver          = 'lsqr'; % 'lsqr' or 'bicg'
par.kernels         = 'FF'; % 'FF' or 'ray' < option to use finite freq kernels %%%%

par.synth_test      = 0;  % this will stop actual inversion

par.crust_corr      = 0; % option to do crustal correction
par.wtdata          = 1; % extra QCs found in wtdata.m function

par.Rdvpdvs         = 0; % 0 or 0.55 % factor for scaling dlnVp to dlnVs (0 means no scaling this way). R = dlnVp/dlnVs, so <1
% damp dVp and dVs smoothnesses together
par.VpVs_smooth_same = 2; % 0 to ignore this, otherwise numebr is weight

par.squeeze         = 0; % 0=no, 1=yes, 0<frac<1=force 2D for first frac of runs
par.zsqz            = 250; %180% squeezing depth in km, if -ive then sqz to below this depth

par.saveopt         = 0;

%% regularisation parms
%{
if par.t_ts == 1
par.damp            = 40;  %3 % scaling for whole damp mat: more damping parms found in dampmat.m  %%%%%
	par.damp_evt    = 0.01; %1 % absolute model damping to use for evt terms
	par.damp_stn    = 15;    %1 % absolute model damping to use for stn terms

par.smooth          = 40; %3 %  scaling for whole smth mat: good result with 1 - %%%%%%
    par.smth_zvh    = 0.5; %0.3 %  Ratio of vertical to horizontal smoothing
elseif par.t_ts == 2
par.damp            = 5e-3;  %3 % scaling for whole damp mat: more damping parms found in dampmat.m
	par.damp_evt    = 0.01; %1 % absolute model damping to use for evt terms
	par.damp_stn    = 1;    %1 % absolute model damping to use for stn terms

par.smooth          = 10.6; %0.1%3 %  scaling for whole smth mat: good result with 1 - 
    par.smth_zvh    = 0.5; %0.3 %  Ratio of vertical to horizontal smoothing
end
%}
if par.t_ts == 1
par.damp            = 2;  %3 % scaling for whole damp mat: more damping parms found in dampmat.m  %%%%%
	par.damp_evt    = 0.01; % 0.01 1 % absolute model damping to use for evt terms
	par.damp_stn    = 1;    %1 1 % absolute model damping to use for stn terms

par.smooth          = 4.5; %3 %  scaling for whole smth mat: good result with 1 - %%%%%%
    par.smth_zvh    = 0.5; %0.3 %  Ratio of vertical to horizontal smoothing
par.dampx = 1;
elseif par.t_ts == 2
par.damp            = 5e-3;  %3 % scaling for whole damp mat: more damping parms found in dampmat.m
	par.damp_evt    = 0.01; %1 % absolute model damping to use for evt terms
	par.damp_stn    = 1;    %1 % absolute model damping to use for stn terms

par.smooth          = 10.6; %0.1%3 %  scaling for whole smth mat: good result with 1 - 
    par.smth_zvh    = 0.5; %0.3 %  Ratio of vertical to horizontal smoothing
end

% more smoothing of certain depths
par.depsmth = [];
%                   minz maxz smooth_wt_multiple
par.depsmth(1,:) = [0    50    2]; % 2x extra damp for shallow nodes
par.depsmth(2,:) = [500  1000  2]; % 2x extra damp for deep nodes

par.scalereg        = 0; %       option to scale regularisation matrices in make_F_f

par.age_smooth      = 0; % option to smooth points together based on age bins - only relevant for oceans - ignore otherwise

% ==== more damping parms in make_dampmat.m ===%

%% plotting parms - zero for speed
par.plot_crustcorr  = 0;
par.plot_data       = 0;
par.plot_inmodel    = 1;
par.plot_synout     = 1;
par.plot_outputs    = 1;
par.plot_Zslices    = 1;
par.plot_hitq       = 0;
par.plot_raypath    = 0;
par.plot_everyNtime = 0; % N>0 will plot each N iterations
 
%% model space parms
par.origin  = [-9 -136];   % model origin [lat,lon]
mstruct = defaultm('mercator');
mstruct.origin = [par.origin 0];
mstruct = defaultm( mstruct ); 
par.map_proj = mstruct;

par.dh      = 30; %80   % horizontal spaces for nodes
par.dh_max  = 30; %100   % horizontal spacing to which it scales, at twice array aperture

par.dz      = 30; %35   % vertical spacing for nodes beneath moho
par.dz_max  = 30; %40   % vertical spacing to which it scales, at model base

par.zmax    = 320; %300 % maximum vertical distance
par.zmin    = 6;  %40  % top of model (above this must be accounted for by crust
 
% kernel integration calculation parms
par.kidd    = 5;      % along-ray segment lengths to calc at
 
%% starting model parms
par.Vavmod  = 1; % starting model: 1=AK135,2=Raj+Ferris,3=STW105
    
%% vertical section parameters
    % give ends of vertical cross sections to compute. As many as you
    % want...
    zsectxy(:,:,1) = [-133.55 -8; -132.45 -4];
    %zsectxy(:,:,1) = [-175 63; -125 63]; % [lon1,lat1; lon2,lat2];
    %zsectxy(:,:,2) = [-155 50; -155 75];

par.zsectxy = zsectxy;

    
%% synth model parms
par.synth_noisy     = 1;  % option to include 'realistic' noise in synthetic data
sym.opt  = 'checker'; % 'checker' or 'custom'

% custom structure
sym.acx = [-580	 -530  -470  -410  -310  -370   130  ]           ;% vector of x-coords for anomaly centres
sym.acy = [  40	  190   325   470  -240  -410   -50  ]           ;% vector of y-coords for anomaly centres
sym.acz = [  70	   70    70    70    70    70   180  ]           ;% vector of z-coords for anomaly centres

sym.awx = [  80	   60    80    80    80    60   100]           ;% vector of x-widths for anomalies (km)
sym.awy = [ 140	  140   120   150   150   140  1000]           ;% vector of y-widths for anomalies (km)
sym.awz = [  90	   90    90    90    90    90   300]           ;% vector of z-widths for anomalies (km)

sym.adval=[  5	   5    5    5    5    5    -2]           ;% vector of pertubations of anomalies (percent)

% checker structure
sym.dval = 5;                        % in percent!
if par.t_ts==2, sym.dval = 300; end

sym.noise = 0.2;                     % 0.05 % std. of gaussian perturbations to data values (s)
sym.estatic_sd = 0.03;               % 0.5 % std. of gaussian perturbations input event terms (s)
sym.sstatic_sd = 0.05;               % 0.5 % std. of gaussian perturbations input station terms (s)

par.sym = sym;

    
%     % if 2D and checker, make 2D checkers...
%     if par.force2D && strcmp(sym.opt,'checker')
% 	sym.acx = [   0    0    0    0    0    0    0    0     0    0    0    0    0    0 ]           ;% vector of x-coords for anomaly centres
%     sym.acy = [-235 -235  -60  -60   60   60  200  200  -130 -130    0    0  120  120 ]           ;% vector of y-coords for anomaly centres
%     sym.acz = [  75  155   75  155   75  155   75  155    75  155   75  155   75  155 ]           ;% vector of z-coords for anomaly centres
%     
%     sym.awx = [ 1e3  1e3  1e3  1e3  1e3  1e3  1e3  1e3   1e3  1e3  1e3  1e3  1e3  1e3 ]           ;% vector of x-widths for anomalies (km)
%     sym.awy = [ 150  150   90   90   90   90  120  120    90   90   90   90   90   90 ]           ;% vector of y-widths for anomalies (km)
%     sym.awz = [  50   50   50   50   50   50   50   50    50   50   50   50   50   50 ]           ;% vector of z-widths for anomalies (km)
% 
%     sym.adv = [  -5    5    5   -5   -5    5    5   -5     0    0    0    0    0    0 ]           ;% vector of v pertubations of anomalies (percent)
%     sym.aa  = [   0    0    0    0    0    0    0    0    -5    5    5   -5   -5    5 ]           ;% vector of anisotropy values of anomalies (percent)
%     sym.opt  = 'custom'; % 'checker' or 'custom'
%     
%     par.sym = sym;
%     end
    
% %     
% %% RUN!
% if ~exist('bokachoda','var')
% run_inversion
% end
% 
% %
%synth_test
