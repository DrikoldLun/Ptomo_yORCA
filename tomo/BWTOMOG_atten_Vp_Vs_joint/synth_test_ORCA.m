function semb_integral = syn_test_ORCA(syn_datafile,damp,smooth)
close all
global F f G res wd dbname
% profile on
%wd = '~/Dropbox/MATLAB/atten_tomog_rays';
try
    shutdown_seizmo;
catch
    1;
end
addpath('matguts','plotting','function','seizmo');


fprintf('\n=========== RUNNING INVERSION ===========\n\n')

%% ===================== NAME OF PROJECT =====================
dbname = 'yORCA';

%% ALL PARMS ESTABLISHED IN PARMFILE
% consider moving frequently changed ones to this script
fprintf('>  Establishing parameters\n')
%cd(wd);
run([dbname,'/PARMS_syn'])
basename = split(datfile,'/');
type = replace(basename{end},'.dat','');
disp(type)
startup_seizmo;
par.damp = damp;
par.smooth = smooth;

if par.t_ts == 1
    Kfile = [dbname,'/K_',par.phasecomp,'_v_',par.kernels,'_',type,'_angle',num2str(par.map_proj.origin(3)),'.mat'];
    modfile = 'model_v';
    parfile = 'par_v';
elseif par.t_ts == 2
    Kfile = [dbname,'/K_',par.phasecomp,'_q_',par.kernels];
    modfile = 'model_q';
    parfile = 'par_q';
end

%% read data

    %datfile = {'yORCA/data/syn_dT_PZ.dat'};
    datfile = {syn_datafile};
    fprintf('>  Reading in data and station details\n')
    [data,par] = read_data(datfile,stafile,crustfile,par);

%% weight data
if par.wtdata
    fprintf('>  Weighting data\n')
    data.ray.wt = wtdata(data);
else
    fprintf('>  No data weighting\n')
    data.ray.wt = ones(size(data.ray.d));
end

%% setup tomo geom
fprintf('>  Setting up geometry\n')
[ data,par ] = setup_geom( data,par );

%% Crustal correction
data = crustal_correction(data,par);

%% ================  setup synthetic model  ================  
synth_model = make_synth_model(par,data);

if par.plot_inmodel
    fprintf('>  Plotting input model\n')
    [plot_simodel] = conv2plotable(synth_model,par);
    plot_Hmaps(plot_simodel,par,2,par.saveopt)
    plot_Zmaps_synth(plot_simodel,par,data,1,par.saveopt)
end

%% =========== Make Kernel =============================================
try
    load(Kfile);
    if length(K.n_indx)~=data.ray.nrays || max(K.n_indx{1}) > par.nmodel
        error('Loaded K not appropriate for this data... build K again')
    end
catch err
    disp(err)
    [K,data] = make_K( par, data );
    save (Kfile,'K','-v7.3');
end

%{
if par.build_K == 0
    load(Kfile);
    if length(K.n_indx)~=data.ray.nrays || max(K.n_indx{1}) > par.nmodel
        error('Loaded K not appropriate for this data... build K again')
    end
elseif par.build_K == 1
    [K,data] = make_K( par, data );
    save (Kfile,'K','-v7.3');
end
%}


%% Hit quality?
par = calc_hitcq( data,par,K );
if par.plot_hitq
    plot_hitq(par);
end

%% Raypaths
if par.plot_raypath
    plot_raypaths(par,data);
end

fprintf('\nSTARTING SYNTHETIC TEST\n')
% profile clear
% profile on

%% Calc synthetic data
fprintf('>  Calculating synthetic data\n')
% data kernel
G = make_G( K,data,par );
% make synthetic data
d_synth = G*synth_model.mval;
% event and station terms
% noise
if par.synth_noisy
    for i = 1:size(d_synth)
        randn('seed',i);
        d_synth(i) = d_synth(i) + random('norm',0,data.ray.sd(i),1);
    end
    %d_synth = d_synth + random('norm',0,par.sym.noise,size(d_synth));
end

[d_synth,~] = add_static_terms(d_synth,G,data,synth_model);

%% Plot data
%{
if par.plot_data
    d_real = data.ray.d; % save the real data
    data.ray.d = d_synth; % put synthetic data into ray struct to plot
    fprintf('>  Plotting synthetic data\n')
    plot_data(data,par,0) % plot synthetic data
    data.ray.d = d_real; % put real data back in
end
%}

%% ==================== DO INVERSION ================== %%
%% ============ copy from here for real thing ========= %%

[model,par] = make_start_model(par,data);

%% data kernel (redundant, but want to have the full inversion structure in here
% G = make_G( K,data,par );
G_0 = G;

tic
count = 0;
solved = false;
while solved == false % loop if you will do squeezing
count = count+1;

%% Predicted data
G = G_0;
d_pred = G*model.mval;
[d_pred,G] = add_static_terms(d_pred,G,data,model);

%% calculate residual and do inverse problem
d_use = d_synth;
res = d_use - d_pred;

%% ADD REGULARISATION
[ F,f ] = make_F_f( G,res,data,par );

%% SOLVE
if strcmp(par.solver,'lsqr')
    fprintf('>  Solving F*dm = f using LSQR\n')
    [dm,flag,relres,iter,resvec] = lsqr( F, f, 1e-5, 6000 );
    if iter==6000, fprintf('Warning LSQR not converging\n'); end
elseif strcmp(par.solver,'bicg')
    fprintf('>  Solving F*dm = f using biconjugate gradient method\n')
    [ dm ] = solve_bicg( F, f, 7e-5, 1e4 );
end

%% UPDATE
if par.PS == 3, sx = 2; else, sx = 1; end % account for extra station terms
[ model ] = model_update( model, dm, par.nmodel,data.evt.nevts,data.stn.nstas*sx);

toc
solved = true;
end

par.si_model = synth_model;
par.so_model = model;
[par.semb,npts] = zelt_semblance_anis(par);

%% Final residual
d_pred = G_0*model.mval;
[d_pred,G] = add_static_terms(d_pred,G_0,data,model);
res = d_use - d_pred;


%% ================== Plotting ==================
plot_sta_static(par,data,model,synth_model)
plot_results_info(par,model,data,d_use,res,synth_model)
% profile viewer
if par.plot_synout
    [plot_somodel] = conv2plotable(model,par);
    plot_Hmaps(plot_somodel,par,3,par.saveopt)
	plot_Zmaps_synth(plot_somodel,par,data,2,par.saveopt)
end

%{
% save figs
resultdir = sprintf('%ssyn/%d/%.1f_%.1f_%d',par.figdir,length(d_synth),par.damp,par.smooth,length(d_synth));
if ~exist(resultdir,'dir')
    mkdir(resultdir)
end
save tmp
figure(87)
print(gcf,'-dpng',fullfile(resultdir,'Hslice'))
figure(14)
print(gcf,'-dpng',fullfile(resultdir,'Vslice1'))
figure(15)
print(gcf,'-dpng',fullfile(resultdir,'Vslice2'))
figure(88)
print(gcf,'-dpng',fullfile(resultdir,'staavg'))
figure(58)
print(gcf,'-dpng',fullfile(resultdir,'error'))
figure(23)
print(gcf,'-dpng',fullfile(resultdir,'data'))
figure(57)
print(gcf,'-dpng',fullfile(resultdir,'stats'))
%}

if strcmp(par.sym.opt,'checker')
    semb = par.semb;
    save('semblance','semb','npts')
end

semb_integral = mean(semb(find(par.hitq>=0.3)));
%save par_synth par

return
end