function [K,data] = make_K( par, data )
fprintf('>  Making K values...  \n')

%close all
tic
% stn, ray, node data
M = length(par.mx);
N = length(data.ray.d);

% frequency-dependent rmax
cfrmax = [ 0.05 0.1 0.3 0.5 1.2 5;... % for P
          0.04 0.05 0.1 0.4 0.6 5];   % for S
rrmax  = [ 220  180 160 130 100 10;...% for P
            250  230 180 120 100 10]; % for S



%------------- main loop through all rays -----------------------

sta_num = data.ray.sta_num;
slats = data.stn.lat;
slons = data.stn.lon;
selvs = data.stn.elv;

bazs = data.ray.baz;
ray_p = data.ray.p;
cf = data.ray.cf;
trdist = data.ray.trdist;
pors = data.ray.ph; 

Si = cell(N,1);
Sj = cell(N,1);
S  = cell(N,1);

if par.plot_raypath
rpath.rx = nan(100,N);
rpath.ry = nan(100,N);
rpath.rz = nan(100,N);
rpath.rt = nan(100,N);
end

pltry = 0;
h = waitbar(0,'Progress through ray tracing');
for iray = 1:N       % or switch to parfor to speed up
        
    % When 3-D load ray here
    
    snum = sta_num(iray); 
    % segment ends
    [rx,ry,rz,dr] = ray_trace_1D(par, ray_p(iray), bazs(iray), slats(snum), slons(snum) ,selvs(snum), false,pors(iray)); % false = no crust
    % segment midpoints
    mx = midpts(rx);
    my = midpts(ry);
    mz = midpts(rz);
    lr = diff(dr);
    
% station-individualised velocity profile
    %     rv = vel_profile(par.PS,rz,slats(snum), slons(snum) ,selvs(snum), false,sseds(snum));
    %     rt = diff(dr)./midpts(rv);
    %     sum(rt)
    mv = vel_profile(pors(iray),mz+selvs(snum));
    rt = lr./mv; % travel time in each segment
    
    % function ray_nodes to find nodes potentially within RF1
    if strcmp(par.kernels,'FF')
	    rayxyzdr = [ rx ry rz dr ];
		rmax = interp1(cfrmax(pors(iray),:),rrmax(pors(iray),:),cf(iray));
        [n_indx, n_vals ] = ray_nodes_ff(rayxyzdr,dr,rmax,par,ray_p(iray),cf(iray),trdist(iray));
        n_vals = -n_vals;
    elseif strcmp(par.kernels,'ray')
        [n_indx, n_lens ] = ray_nodes_box(mx,my,mz,lr,par);
        % KERNEL VALUE IS THE MEAN TIME PER NODE - model is dSLOWNESS
        % + need negative, as positive delay means negative diffTT ????????
        n_vals = n_lens./par.mvav(n_indx);
        
    end
    
    if par.t_ts == 1 % if velocity inversion        
%         if abs(sum(n_vals))<(par.max_z*par.PS/8)
%             pause
%         end
    elseif par.t_ts == 2 % if q inversion
        % KERNEL VALUE IS THE MEAN TSTAR PER NODE - model is dATTENUATION
        % + need negative, as positive delay means negative diff TT ????????
        % use indices for p or s, as relevant
        if par.PS == 3 && par.Rdvpdvs == 0 && pors(iray) == 2
            pix = par.nmodel;
        else
            pix = 0;
        end
        
        mQavi = par.mQav(n_indx + pix);
        n_vals = n_vals./mQavi(:);     
        
%         if sum(n_vals)<0.1
%             pause
%         end
    end

    
    % %----------------- comment out if parallelised -----------------
    %% PLOT RAYS
     if par.plot_raypath
       rpath.rx(1:length(rx),iray) = rx;
       rpath.ry(1:length(ry),iray) = ry;
       rpath.rz(1:length(rz),iray) = rz;
       rpath.rt(1:length(rt),iray) = rt;
       
         if pltry    
         figure(1); hold on
         scatter3(par.mx(n_indx),par.my(n_indx),par.mz(n_indx),500*(0.000001 + n_vals),'r');
         scatter3(par.mx(n_indx),par.my(n_indx),par.mz(n_indx),500*(0.000001 + par.mvav(n_indx)),'r');
         plot3(rx,ry,rz,'-o','LineWidth',2,'color',colour_get(data.ray.d(iray),-3,3))
         set(gca,'zdir','reverse','xlim',[par.min_x par.max_x],'ylim',[par.min_y par.max_y],'zlim',[0 par.max_z])
         title(sprintf('orid %u, sta %u',data.ray.orid(iray),snum))
         
         
         cont = input('Press any key other than "x" to plot next ray. ','s');
         if strcmp(cont,'x'); pltry = 0; end
         scatter3(par.mx(n_indx),par.my(n_indx),par.mz(n_indx),500*(0.000001 + n_vals),'b');
         end
     end
% 
     if mod(iray,50)==1
     waitbar(iray/N,h)
     end
    % %----------------- comment out if parallelised -----------------
  
    Si{iray} = ones(size(n_indx))*iray;
    Sj{iray} = n_indx;
    S{iray}  = n_vals;

end


if par.plot_raypath
    L = max(sum(~isnan(rpath.rx)));
    rpath.rx(L+1:end,:)=[];
    rpath.ry(L+1:end,:)=[];
    rpath.rz(L+1:end,:)=[];
    rpath.rt(L+1:end,:)  =[]; %modify
    save ('rpath','rpath','-v7.3');
end

K.n_indx = Sj;
K.n_vals = S;
delete(h)
pause(0.01)

toc

