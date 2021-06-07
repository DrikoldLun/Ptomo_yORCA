clear all
close all
global F f G res wd dbname
% profile on
%wd = '~/Dropbox/MATLAB/atten_tomog_rays';
try
    shutdown_seizmo;
catch
    1;
end
addpath('matguts','plotting','function','seizmo','Tinycodes');

fprintf('\n=========== RUNNING INVERSION ===========\n\n')

%% ===================== NAME OF PROJECT =====================
dbname = 'yORCA';

%% ALL PARMS ESTABLISHED IN PARMFILE
% consider moving frequently changed ones to this script
fprintf('>  Establishing parameters\n')
%cd(wd);
run([dbname,'/PARMS'])
basename = split(datfile,'/');
type = replace(basename{end},'.dat','');
disp(type)
startup_seizmo;
is2D = true;
typeadd = '';
if is2D
    par.map_proj.origin(3) = -15;
    par.dampx = 1000;
    typeadd = '_2D';
end

if par.t_ts == 1
    Kfile = [dbname,'/K_',par.phasecomp,'_v_',par.kernels,'_',type,'_angle',num2str(par.map_proj.origin(3)),'.mat'];
    modfile = 'model_v';
    parfile = 'par_v';
elseif par.t_ts == 2
    Kfile = [dbname,'/K_',par.phasecomp,'_q_',par.kernels,'.mat'];
    modfile = 'model_q';
    parfile = 'par_q';
end

%% read data
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

%% Plot data
if par.plot_data
    fprintf('>  Plotting data\n')
    plot_data(data,par,1)
end


%% Make starting model - required to put starting model into parm
[model,par] = make_start_model(par,data);

model_1 = model;

if par.plot_inmodel && ~par.synth_test
    [plot_inmodel] = conv2plotable(model,par);
    %plot_Hmaps(plot_inmodel,par,2,par.saveopt)
    %plot_Zmaps(plot_inmodel,par,2,par.saveopt)
end

%% =========== Make Kernel =============================================
try
    load(Kfile);
    if length(K.n_indx)~=data.ray.nrays || max(K.n_indx{1}) > par.nmodel
        %% 
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

%% Synthetic test
if par.synth_test
    synth_test
    return
end


%% ==================== DO INVERSION ================= %%
G = make_G( K,data,par );

%% predicted data from starting model
d_pred = G*model.mval;
[d_pred,G] = add_static_terms(d_pred,G,data,model);

%% calculate residual and do inverse problem
d_use = data.ray.d ;%- data.ray.ccorr;
res = d_use - d_pred;

%% ADD REGULARISATION
[ F,f ] = make_F_f( G,res,data,par );

%% SOLVE
if strcmp(par.solver,'lsqr')
    fprintf('>  Solving F*dm = f using LSQR\n')
    [dm,flag,relres,iter,resvec] = lsqr( F, f, 1e-5, 500 );
    if iter==500, fprintf('Warning LSQR not converging\n'); end
elseif strcmp(par.solver,'bicg')
    fprintf('>  Solving F*dm = f using biconjugate gradient method\n')
    [ dm ] = solve_bicg( F, f, 1e-4, 1e4 );
end

%% UPDATE
[ model ] = model_update( model, dm, par.nmodel,data.evt.nevts,data.stn.nstas);

%% Final residual
d_pred = G*[model.mval;model.estatic;model.sstatic];
res = d_use - d_pred;
resid = norm(res)

fprintf('>  Results summary:\n')
% profile viewer


%% RESULT STATS
vr  = variance_reduction(d_use,d_pred);
wvr = variance_reduction(data.ray.wt.*d_use,data.ray.wt.*d_pred);

fprintf('Final RMS error: %.2f\n',rms(res));
fprintf('Variance reduction = %.2f %%\n',vr);
fprintf('Weighted variance reduction = %.2f %%\n',wvr);

fprintf('RMS event static = %.2f s \n',rms(model.estatic))
fprintf('RMS station static = %.2f s \n',rms(model.sstatic))

%% high hitq res
% [model_hq ] = hiQmodel( model,model_1,par,0.3 );
% d_pred_hq = G*[model_hq.mval;model_hq.estatic;model_hq.sstatic];
% res_hq     = d_use - d_pred_hq;
% vr_hq = variance_reduction(d_use,d_pred_hq);
%
% fprintf('HI-Q Variance reduction = %.2f %%\n',vr_hq);

% %% Plot data and residual
% if par.plot_data
%     d_real = data.ray.d; % save the real data
%     data.ray.d = res; % put synthetic data into ray struct to plot
%     fprintf('>  Plotting residual\n')
%     plot_data(data,par,0) % plot synthetic data
%     data.ray.d = d_real; % put real data back in
%     clone_figure(32,33);
%     plot_data(data,par,0) % plot synthetic data
%
% end


%% plots
if par.plot_outputs
    plot_results_info(par,model,data,d_use,res,model_1)
    plot_model = conv2plotable(model,par);
    plot_Hmaps(plot_model,par,1,par.saveopt)
    if par.plot_Zslices
        plot_Zmaps(plot_model,par,data,par.saveopt)
    end
end
if par.savemodel
    save(['model/model_',type,typeadd,'.mat'],'model','par','data','plot_model','vr','wvr','resid')
    %figure(15)
    %print('-dpng',['/media/lun/easystore/tmp/poster/model/flatness/',type,typeadd,'_Zmap.png'])
    %close(15)
    %figure(67)
    %print('-dpng',['/media/lun/easystore/tmp/poster/model/flatness/',type,typeadd,'_Hslice.png'])
    %close(67)
end
return

%% save results

if par.saveopt
    save([par.tomoresdir,modfile],'model','-v7.3');
    save([par.tomoresdir,parfile],'par')
end

% save figs
resultdir = sprintf('%sobs/%d/%.1f_%.1f_%d',par.figdir,length(d_use),par.damp,par.smooth,length(d_use));
if ~exist(resultdir,'dir')
    mkdir(resultdir)
end
figure(67)
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

% some squeezing test stats...
fprintf('\nSqueezing analysis:\n')
par.niter = 10;
if par.niter>0
maxitersq = ceil(par.squeeze*par.niter);
fprintf('With squeezing to %.0f km, achieve\n', par.zsqz)
fprintf('  %.2f %% of final total var reduction\n',100*vr.all(maxitersq)/vr.all(end))
fprintf('  %.2f %% of final tdiff var reduction\n',100*vr.tdiff(maxitersq)/vr.tdiff(end))
fprintf('  %.2f %% of final dT var reduction\n',100*vr.dT(maxitersq)/vr.dT(end))
fprintf('  %.2f %% of final wtvar reduction\n',100*wvr.all(maxitersq)/wvr.all(end))
end

indsqz = par.mz>par.zsqz;
fprintf('Norm of dV structure below zsqz = %.2f ',norm(model.dv(indsqz)));
fprintf(' = %.2f%% of a total norm of %.2f\n',100*(norm(model.dv(indsqz))./norm(model.dv))^2,norm(model.dv));
fprintf('Norm of A structure below zsqz = %.2f ',norm(model.aa(indsqz)));
fprintf(' = %.2f%% of a total norm of %.2f\n',100*(norm(model.aa(indsqz))./norm(model.aa))^2,norm(model.aa));

%% plots etc.
load('semblance');
par.semb = semb;
