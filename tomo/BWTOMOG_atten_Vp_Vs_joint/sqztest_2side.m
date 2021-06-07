function sqztest_2side(syn_datafile,is2D)
%clear all
close all
global G wd dbname
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
%{
if is2D
    par.dampx = 1000;
    par.map_proj.origin(3) = -15;
end
%}
basename = split(syn_datafile,'/');
type = replace(basename{end},'.dat','');
disp(type)
load(['model/model_',type,'.mat'])
static = [model.estatic;model.sstatic];
clear model; clear data; clear par; clear plot_model; clear resid; clear vr; clear wvr;
run([dbname,'/PARMS'])
datfile = {syn_datafile};
direction = 'both';
startup_seizmo;

%% read data
fprintf('>  Reading in data and station details\n')
[data,par] = read_data(datfile,stafile,crustfile,par);

if par.t_ts == 1
    Kfile = [dbname,'/K_',par.phasecomp,'_v_',par.kernels,'_',type,'_angle',num2str(par.map_proj.origin(3)),'.mat'];
    modfile = 'model_v';
    parfile = 'par_v';
elseif par.t_ts == 2
    Kfile = [dbname,'/K_',par.phasecomp,'_q_',par.kernels];
    modfile = 'model_q';
    parfile = 'par_q';
end

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

%% Make starting model - required to put starting model into parm
[model,par] = make_start_model(par,data);

model_1 = model;

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

%% Hit quality?
par = calc_hitcq( data,par,K );


%% ==================== DO INVERSION ================= %%
G = make_G( K,data,par );

%% predicted data from starting model
d_pred = G*model.mval;
[d_pred,G] = add_static_terms(d_pred,G,data,model);
nrays = data.ray.nrays;
WG = spdiags(data.ray.wt,0,nrays,nrays)*G;

%% calculate residual and do inverse problem
d_use = data.ray.d ;%- data.ray.ccorr;
res = d_use - d_pred;
res0 = res;

%% ADD REGULARISATION
[ F,f ] = make_F_f( G,res,data,par );
num = length(par.zz);

is = 1;
for i = 1:num
zs = par.zz(i);
for j = i:num
zd = par.zz(j);
inds = find(par.mz>=zs&par.mz<=zd);
ns = length(inds);

[ Fs,fs] = make_F_f_dep( G,res0,data,par,[zs zd],direction );
colind = [inds',par.nmodel+1:size(Fs,2)];
Fs = Fs(:,colind);
[ H,h ] = make_H_h(Fs,fs,length(inds),static);

%{
A = Fs'*Fs;
WG = spdiags(data.ray.wt,0,nrays,nrays)*G;
WG = WG(:,colind);
Gg = A\WG';
R = Gg*WG;
reso = full(diag(R));
squeezetest.freedom(is) = sum(reso(1:ns));
squeezetest.reso{is} = reso(1:ns);
fprintf('Reso: %.1f\n',sum(reso(1:ns)));
%}

%% m0 solver
if strcmp(par.solver,'lsqr')
    fprintf('>  Solving F*dm = f using LSQR\n')
    [dm0,flag,relres,iter,resvec] = lsqr( F, f, 1e-5, 1000 );
    if iter==500, fprintf('Warning LSQR not converging\n'); end
elseif strcmp(par.solver,'bicg')
    fprintf('>  Solving F*dm = f using biconjugate gradient method\n')
    [ dm0 ] = solve_bicg( F, f, 1e-4, 1e4 );
end
[ model0 ] = model_update( model, dm0, par.nmodel,data.evt.nevts,data.stn.nstas);

%% m1 solver
%{
if strcmp(par.solver,'lsqr')
    fprintf('>  Solving F*dm = f using LSQR\n')
    [dm1,flag,relres,iter,resvec] = lsqr( Fs, fs, 1e-5, 1000 );
    if iter==500, fprintf('Warning LSQR not converging\n'); end
elseif strcmp(par.solver,'bicg')
    fprintf('>  Solving F*dm = f using biconjugate gradient method\n')
    [ dm1 ] = solve_bicg( Fs, fs, 1e-4, 1e4 );
end
%}
dm1 = H\h;
dm1 = [zeros(sum(par.mz<zs),1);dm1(1:ns);zeros(sum(par.mz>zd),1);dm1(length(inds)+1:length(inds)+length(static))];
[ model1 ] = model_update( model, dm1, par.nmodel,data.evt.nevts,data.stn.nstas);

%% Final residual
d_pred = G*[model1.mval;model1.estatic;model1.sstatic];
res = d_use - d_pred;

%% m2 solver
F2 = F; f2 = [res;f(length(res)+1:end)];
if strcmp(par.solver,'lsqr')
    fprintf('>  Solving F*dm = f using LSQR\n')
    [dm2,flag,relres,iter,resvec] = lsqr( F2, f2, 1e-5, 1000  );
    if iter==500, fprintf('Warning LSQR not converging\n'); end
elseif strcmp(par.solver,'bicg')
    fprintf('>  Solving F*dm = f using biconjugate gradient method\n')
    [ dm2 ] = solve_bicg( F2, f2, 1e-4, 1e4 );
end
[ model2 ] = model_update( model, dm2, par.nmodel,data.evt.nevts,data.stn.nstas);

fprintf('>  Results summary:\n')
% profile viewer


%% RESULT STATS
fprintf('squeezing depth: %.1fkm-%.1fkm\n',zs,zd);
vr  = variance_reduction(d_use,d_pred);
wvr = variance_reduction(data.ray.wt.*d_use,data.ray.wt.*d_pred);

fprintf('Final RMS error: %.2f\n',rms(res));
fprintf('Variance reduction = %.2f %%\n',vr);
fprintf('Weighted variance reduction = %.2f %%\n',wvr);

fprintf('RMS event static = %.2f s \n',rms(model1.estatic))
fprintf('RMS station static = %.2f s \n',rms(model1.sstatic))

%% save results
%{
if par.saveopt
    save([par.tomoresdir,modfile],'model','-v7.3');
    save([par.tomoresdir,parfile],'par')
end
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
%}

%indsqz = par.mz>par.zsqz;
%fprintf('Norm of dV structure below zsqz = %.3f ',norm(model1.mval(indsqz)));
%fprintf(' = %.3f%% of a total norm of %.3f\n',100*(norm(model1.mval(indsqz))./norm(model1.mval)),norm(model1.mval));
%fprintf('m2 norm reduction = %.3f%%\n',100*(1-norm(model2.mval)/norm(model0.mval)));
%fprintf('Norm of A structure below zsqz = %.2f ',norm(model.aa(indsqz)));
%fprintf(' = %.2f%% of a total norm of %.2f\n',100*(norm(model.aa(indsqz))./norm(model.aa))^2,norm(model.aa));

%% plots

if par.plot_outputs
    %plot_results_info(par,model,data,d_use,res,model_1)
    plot_model = conv2plotable(model1,par);
    plot_Hmaps(plot_model,par,1,par.saveopt)
    figure(67)
    print('-dpng',['/media/lun/easystore/tmp/IM4.27/sqztest/2side/',num2str(zs),'_',num2str(zd),'_Hslice.png'])
    close(67)
    if par.plot_Zslices
        plot_Zmaps(plot_model,par,data,par.saveopt)
        figure(15)
        print('-dpng',['/media/lun/easystore/tmp/IM4.27/sqztest/2side/',num2str(zs),'_',num2str(zd),'_Zmap.png'])
        close(15)
        close(14)
    end
end

squeezetest.model{is} = model1;
squeezetest.is(is) = i;
squeezetest.id(is) = j;
squeezetest.nmodel(is) = ns;
squeezetest.zs(is) = zs;
squeezetest.zd(is) = zd;
squeezetest.vr(is) = vr;
squeezetest.wvr(is) = wvr;
squeezetest.resid(is) = norm(res);
squeezetest.m2norm(is) = norm(model2.mval);
squeezetest.m2reduction(is) = 100*(1-norm(model2.mval)/norm(model0.mval));
is = is + 1;
end
end
if ~is2D
    ofile = ['sqztest_2side/sqz2side_',type,'.mat'];
else
    ofile = ['sqztest_2side/sqz2side_',type,'_2D.mat'];
end
save(ofile,'squeezetest','par')

end
