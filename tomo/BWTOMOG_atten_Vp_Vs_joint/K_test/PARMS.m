clear par 
% project details
global dbdir
dbdir = 'K_test/'; % include final slash

if ~exist('par','var')==1, par=struct([]); end

par(1).phasecomp = 'PZ'; % 'PZ' or 'ST'
par.dattype = 'dT'; % 'dT' or 'dtstar'

%% ---------------------------------
par.PS = find(strcmp({'P','S'},par.phasecomp(1))); % 1 for P, 2 for S
par.t_ts = find(strcmp({'dT','dtstar'},par.dattype)); % 1 for dT, 2 for dtstar

%% results directories
% run([dbdir,dbname,'_startup.m']);
par.figdir = [dbdir,'figs/'];
par.tomoresdir = [dbdir,'results/'];
par.tomodatdir = [dbdir,'data/'];

%datfile = sprintf('%sdata_%s_%s.dat',par.tomodatdir,par.dattype,par.phasecomp);
%stafile  = sprintf('%sstations_%s.dat',par.tomodatdir,par.phasecomp);
[datfile,stafile,par.phases] = parse_phasecomp(par.tomodatdir,par.phasecomp,par.dattype);
crustfile  = sprintf('%sstations_crust.dat',par.tomodatdir);

%% inversion running parms - zero for speed
REBUILD             = 1;
par.redo_crust_corr = REBUILD; % <== REDO IF DATA/MODEL CHANGES
par.build_K         = REBUILD; % <== REDO IF DATA/MODEL CHANGES
par.build_smooth    = REBUILD; % <== REDO IF DATA/MODEL CHANGES
par.solver          = 'lsqr'; % 'lsqr' or 'bicg'
par.kernels         = 'ray'; % 'FF' or 'ray' < option to use finite freq kernels

par.synth_test      = 1;  % this will stop actual inversion

par.crust_corr      = 0; % option to do crustal correction
par.wtdata          = 1; % extra QCs found in wtdata.m function

par.squeeze         = 0; % 0=no, 1=yes, 0<frac<1=force 2D for first frac of runs
par.zsqz            = 175; %180% squeezing depth in km, if -ive then sqz to below this depth

par.saveopt         = 0;

%% regularisation parms
par.age_smooth      = 0; % option to smooth points together based on age bins

if par.t_ts == 1
par.damp            = 10;  %3 % scaling for whole damp mat: more damping parms found in dampmat.m
	par.damp_evt    = 0.01; %1 % absolute model damping to use for evt terms
	par.damp_stn    = 3;    %1 % absolute model damping to use for stn terms

par.smooth          = 12; %3 %  scaling for whole smth mat: good result with 1 - 
    par.smth_zvh    = 0.5; %0.3 %  Ratio of vertical to horizontal smoothing
elseif par.t_ts == 2
par.damp            = 5e-3;  %3 % scaling for whole damp mat: more damping parms found in dampmat.m
	par.damp_evt    = 0.01; %1 % absolute model damping to use for evt terms
	par.damp_stn    = 3;    %1 % absolute model damping to use for stn terms

par.smooth          = .1; %3 %  scaling for whole smth mat: good result with 1 - 
    par.smth_zvh    = 0.5; %0.3 %  Ratio of vertical to horizontal smoothing
end

par.scalereg        = 0; %       option to scale regularisation matrices in make_F_f
% ==== more damping parms in make_dampmat.m ===%

%% plotting parms - zero for speed
par.plot_crustcorr  = 0;
par.plot_data       = 0;
par.plot_inmodel    = 0;
par.plot_synout     = 1;
par.plot_outputs    = 1;
par.plot_Zslices    = 1;
par.plot_hitq       = 0;
par.plot_raypath    = 0;
par.plot_everyNtime = 0; % N>0 will plot each N iterations
 
%% model space parms
    par.origin  = [0.1 0.1];   % model origin [lat,lon]
    mstruct = defaultm('mercator');
    mstruct.origin = [par.origin 0];
    mstruct = defaultm( mstruct ); 
    par.map_proj = mstruct;
 
    par.dh      = 4; %30   % horizontal spaces for nodes
    par.dh_max  = 4; %80   % horizontal spacing to which it scales, at twice array aperture

    par.dz      = 4; %35   % vertical spacing for nodes beneath moho
    par.dz_max  = 4; %40   % vertical spacing to which it scales, at model base
 
    par.zmax    = 200; %300 % maximum vertical distance
    par.zmin    = 4;  %40  % top of model (above this must be accounted for by crust
 
% kernel integration calculation parms
    par.kidd    = 2;      % along-ray segment lengths to calc at
 
%% starting model parms
    par.Vavmod  = 2; % starting model: 1=AK135,2=Raj+Ferris,3=STW105
    
%% synth model parms
    par.synth_noisy     = 1;  % option to include 'realistic' noise in synthetic data
    sym.opt  = 'halves'; % 'checker' or 'custom'
    
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
        
    sym.noise = 0.0001;                     % 0.05 % std. of gaussian perturbations to data values (s)
    sym.estatic_sd = 0.00;               % 0.5 % std. of gaussian perturbations input event terms (s)
    sym.sstatic_sd = 0.00;               % 0.5 % std. of gaussian perturbations input station terms (s)
    
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
