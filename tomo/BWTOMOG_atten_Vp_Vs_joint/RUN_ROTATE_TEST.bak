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
addpath('matguts','plotting','function','seizmo');

fprintf('\n=========== RUNNING INVERSION ===========\n\n')

%% ===================== NAME OF PROJECT =====================
dbname = 'yORCA';

%% ALL PARMS ESTABLISHED IN PARMFILE
% consider moving frequently changed ones to this script
fprintf('>  Establishing parameters\n')
%cd(wd);
run([dbname,'/PARMS'])
datfile = 'yORCA/data/PZmix.dat'; %2,5
%datfile = 'yORCA/data/Z_mag_5.5_10.0.dat'; %1,5
%datfile = 'yORCA/data/PZ_direct_num.dat';%2,5
%datfile = 'yORCA/data/P_mag_6.0_10.0.dat'; %2,5
basename = split(datfile,'/');
type = replace(basename{end},'.dat','');
disp(type)
K_dir = ['rotate_K_',type];
if ~isdir(K_dir)
    mkdir(K_dir)
end
startup_seizmo;
%rotate_test.angle = -[0:5:180]';
rotate_test.angle = -[0:5:30]';
num = length(rotate_test.angle);
rotate_test.resid = zeros(num,1);
rotate_test.vr = zeros(num,1);
rotate_test.wvr = zeros(num,1);

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

for i = 1:num
par.map_proj.origin(3) = rotate_test.angle(i);
if par.t_ts == 1
    Kfile = [K_dir,'/K_',par.phasecomp,'_v_',par.kernels,'_',type,'_angle',num2str(par.map_proj.origin(3)),'.mat'];
    modfile = 'model_v';
    parfile = 'par_v';
elseif par.t_ts == 2
    Kfile = ['rotate_K','/K_',par.phasecomp,'_q_',par.kernels,'.mat'];
    modfile = 'model_q';
    parfile = 'par_q';
end

%% setup tomo geom
fprintf('>  Setting up geometry\n')
[ data,par ] = setup_geom( data,par );


%% Crustal correction
data = crustal_correction(data,par);

%% Make starting model - required to put starting model into parm
[model,par] = make_start_model(par,data);

model_1 = model;

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

rotate_test.resid(i) = norm(res);
rotate_test.vr(i) = vr;
rotate_test.wvr(i) = wvr;

end

ofile = ['rotate_test_',type,'.mat'];
save(ofile,'rotate_test')

%{
%% plots
if par.plot_outputs
    plot_results_info(par,model,data,d_use,res,model_1)
    plot_model = conv2plotable(model,par);
    plot_Hmaps(plot_model,par,1,par.saveopt)
    if par.plot_Zslices
        plot_Zmaps(plot_model,par,data,par.saveopt)
    end
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
%}
