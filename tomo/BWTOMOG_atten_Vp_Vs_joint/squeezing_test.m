function squeezing_test(syn_datafile,direction,is2D)
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
%typeadd = '';
%{
if is2D
    par.dampx = 1000;
    par.map_proj.origin(3) = -15;
    typeadd = '_2D';
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
load(['sqztest_',type,'_',direction,'.mat'])
squeezetestbak = squeezetest;
clear squeezetest

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
%{
squeezetest.zsqz = zeros(num,1);
squeezetest.vr = zeros(num,1);
squeezetest.wvr = zeros(num,1);
squeezetest.resid = zeros(num,1);
squeezetest.m2norm = zeros(num,1);
squeezetest.m2reduction = zeros(num,1);
%}

if strcmp(direction,'down')
indep = length(par.zz):-1:1;
elseif strcmp(direction,'up')
indep = 1:length(par.zz);
end

is = 1;
for i = indep
par.zsqz = par.zz(i);
if strcmp(direction,'down')
ns = sum(par.mz<=par.zsqz);
inds = find(par.mz<=par.zsqz);
elseif strcmp(direction,'up')
ns = sum(par.mz>=par.zsqz);
inds = find(par.mz>=par.zsqz);
end
%{
constant = 1000;
S = sparse(ns,size(F,2));
s = zeros(ns,1);
for is = 1:ns
    S(is,inds(is)) = constant;
end
%}
[ Fs,fs] = make_F_f_dep( G,res0,data,par,par.zsqz,direction );
colind = [inds',par.nmodel+1:size(Fs,2)];
Fs = Fs(:,colind);
[ H,h ] = make_H_h(Fs,fs,length(inds),static);
%{
Fs = [Fs;S];
fs = [fs;s];
%}
%{
A = Fs'*Fs;
WG = spdiags(data.ray.wt,0,nrays,nrays)*G;
WG = WG(:,colind);
Gg = A\WG';
R = Gg*WG;
reso = full(diag(R));
squeezetest.freedom(i) = sum(reso(1:length(inds)));
squeezetest.reso{i} = reso(1:length(inds));
squeezetest.nmodel(i) = length(inds);
fprintf('Reso: %.1f\n',sum(reso(1:length(inds))));
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

%{
%% m1 solver
if strcmp(par.solver,'lsqr')
    fprintf('>  Solving F*dm = f using LSQR\n')
    [dm1,flag,relres,iter,resvec] = lsqr( H, h, 1e-5, 1000 );
    if iter==500, fprintf('Warning LSQR not converging\n'); end
elseif strcmp(par.solver,'bicg')
    fprintf('>  Solving F*dm = f using biconjugate gradient method\n')
    [ dm1 ] = solve_bicg( H, h, 1e-4, 1e4 );
end
%}
dm1 = H\h;
if strcmp(direction,'down')
    dm1 = [dm1(1:length(inds));zeros(par.nmodel-length(inds),1);dm1(length(inds)+1:length(inds)+length(static))];
elseif strcmp(direction,'up')
    dm1 = [zeros(par.nmodel-length(inds),1);dm1(1:length(inds)+length(static))];
end
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
fprintf('squeezing depth: %.1fkm\n',par.zsqz);
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

indsqz = par.mz>par.zsqz;
fprintf('Norm of dV structure below zsqz = %.3f ',norm(model1.mval(indsqz)));
fprintf(' = %.3f%% of a total norm of %.3f\n',100*(norm(model1.mval(indsqz))./norm(model1.mval)),norm(model1.mval));
fprintf('m2 norm reduction = %.3f%%\n',100*(1-norm(model2.mval)/norm(model0.mval)));
%fprintf('Norm of A structure below zsqz = %.2f ',norm(model.aa(indsqz)));
%fprintf(' = %.2f%% of a total norm of %.2f\n',100*(norm(model.aa(indsqz))./norm(model.aa))^2,norm(model.aa));

%% plots

if par.plot_outputs
    %plot_results_info(par,model,data,d_use,res,model_1)
    plot_model = conv2plotable(model1,par);
    plot_Hmaps(plot_model,par,1,par.saveopt)
    figure(67)
    print('-dpng',['/media/lun/easystore/tmp/IM4.27/sqztest/',direction,'/',num2str(par.zsqz),'km_',direction,'_Hslice.png'])
    close(67)
    if par.plot_Zslices
        plot_Zmaps(plot_model,par,data,par.saveopt)
        figure(15)
        print('-dpng',['/media/lun/easystore/tmp/IM4.27/sqztest/',direction,'/',num2str(par.zsqz),'km_',direction,'_Zmap.png'])
        close(15)
        close(14)
    end
end

%squeezetest.nmodel(i) = ns;
squeezetest.model{is} = model1;
squeezetest.zsqz(is) = par.zsqz;
squeezetest.vr(is) = vr;
squeezetest.wvr(is) = wvr;
squeezetest.resid(is) = norm(res);
squeezetest.m2norm(is) = norm(model2.mval);
squeezetest.m2reduction(is) = 100*(1-norm(model2.mval)/norm(model0.mval));
%squeezetest.estatic{is} = model1.estatic;
%squeezetest.sstatic{is} = model1.sstatic;
is = is + 1;

figure(1)
clf;
set(gcf,'position',[500 500 500 300])
[AX,H1,H2] = plotyy(squeezetestbak.zsqz,squeezetestbak.wvr,squeezetestbak.zsqz,squeezetestbak.m2norm);
%set([H1;H2],'Marker','.','MarkerSize',24)
set(AX,'position',[0.15,0.1,0.72,0.8])
set(H1,'LineWidth',2);
set(H2,'LineWidth',2);
set(AX,'Xlim',[min(par.zz) max(par.zz)],'FontSize',8,'Xtick',par.zz)


%set(AX(1),'ytick',ceil(min(squeezetest.wvr./squeezetest.nmodel')*500)/500:floor((max(squeezetest.wvr./squeezetest.nmodel')-min(squeezetest.wvr./squeezetest.nmodel'))*500)/2500:max(squeezetest.wvr./squeezetest.nmodel'))
%set(AX(1),'ytick',ceil(min(squeezetest.wvr)):floor((max(squeezetest.wvr)-min(squeezetest.wvr))/5):max(squeezetest.wvr))
axes(AX(1))
hold on
plot(par.zsqz,wvr,'Marker','.','MarkerSize',30,'color','blue')
set(AX(1),'ytick',linspace(65,90,5))
ylabel(AX(1),'$Variance\ Reduction(\%)$','interpreter','latex','Fontsize',12,'FontWeight','bold')
ylim(AX(1),[65 90])

axes(AX(2))
hold on
plot(par.zsqz,norm(model2.mval),'Marker','.','MarkerSize',30,'color',[0.8500, 0.3250, 0.0980])
%set(AX(2),'ytick',ceil(min(squeezetest.m2norm)*500)/500:floor((max(squeezetest.m2norm)-min(squeezetest.m2norm))*500)/2500:max(squeezetest.m2norm))
ylabel(AX(2),'$m2\ norm$','interpreter','latex','Fontsize',12,'FontWeight','bold')
ylim(AX(2),[0.2 0.35])
set(AX(2),'ytick',linspace(0.2,0.35,5))
xlabel('$zsqz(km)$','interpreter','latex','Fontsize',12,'FontWeight','bold')
grid()
%title(titletxt,'interpreter','latex','Fontsize',12,'FontWeight','bold')
figure(1)
print('-dpng',['/media/lun/easystore/tmp/IM4.27/sqztest/',direction,'/',num2str(par.zsqz),'km_',direction,'_sqzcurve.png'])
close(1)

end
if ~is2D
    ofile = ['sqztest/sqztest_',type,'_',direction,'.mat'];
else
    ofile = ['sqztest/sqztest_',type,'_',direction,'_2D.mat'];
end
save(ofile,'squeezetest','par')

end
