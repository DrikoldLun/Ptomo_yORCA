function [ data,par ] = setup_geom( data,par )
% [ data,par ] = setup_geom( data,par )
% setup the geometry of the model space for the inverse problem.



%% ========== define model space ==========================================
[stax,stay] = project_xy( par, data.stn.lat, data.stn.lon);
[stax1,stay1] = project_xy( par, data.stn.alllat, data.stn.alllon);
data.stn.stax = stax;
data.stn.stay = stay;

% plot limits
Wmost_sta = min(data.stn.alllon);
Emost_sta = max(data.stn.alllon);
Nmost_sta = max(data.stn.alllat);
Smost_sta = min(data.stn.alllat);
%Wmost_sta = -135.5;
%Emost_sta = -130;
%Nmost_sta = -4;
%Smost_sta = -8.5;
par.plot_lonlims = round_level([Wmost_sta Emost_sta] + diff([Wmost_sta Emost_sta])*0.1*[-1 1],0.5);
par.plot_latlims = round_level([Smost_sta Nmost_sta] + diff([Smost_sta Nmost_sta])*0.1*[-1 1],0.5);
%par.plot_lonlims = [-136 -130];
%par.plot_latlims = [-8.5 -3.5];

% increase node spacing from centre out and down
maxx = round_level(max(stax1),par.dh); 
minx = round_level(min(stax1),par.dh);
maxy = round_level(max(stay1),par.dh); 
miny = round_level(min(stay1),par.dh);
% lx = maxx-minx;
% ly = maxy-miny;
xx = [minx:par.dh:maxx];
yy = [miny:par.dh:maxy];

% endy = cumsum(round_level(linspace(par.dh,par.dh_max,ceil(lx/(mean([par.dh,par.dh_max])))),5));
endxy = cumsum(round_level(linspace(par.dh,par.dh_max,ceil(par.zmax/(mean([par.dh,par.dh_max])))),5));

% define nodes
xx = [minx-fliplr(endxy) xx maxx+endxy]';
yy = [miny-fliplr(endxy) yy maxy+endxy]';
%zz = [0 par.zmin par.zmin+cumsum(round_level(linspace(par.dz,par.dz_max,floor((par.zmax-par.zmin)/(mean([par.dz,par.dz_max])))),5))]';
zz = [6,40,80,110,140,180,220,260,300]';
% put nodes into vectors
save tmp;
[X,Y,Z] = meshgrid(xx,yy,zz);
    % cycle over y then x then z.
mx = reshape(X,numel(X),1);
my = reshape(Y,numel(X),1);
mz = reshape(Z,numel(X),1);
[mlt, mln] = project_xy(par, mx, my, 'inverse');

dxx = 0.5*([0;diff(xx)]+[diff(xx);0]);
dyy = 0.5*([0;diff(yy)]+[diff(yy);0]);
% dzz = 0.5*([0;diff(zz)]+[diff(zz);0])
dzz = [0;diff(zz)];
[X,Y,Z] = meshgrid(dxx,dyy,dzz);
mdx = reshape(X,numel(X),1);
mdy = reshape(Y,numel(X),1);
mdz = reshape(Z,numel(X),1);

%% save into par structure
par.mx = mx;
par.my = my;
par.mz = mz;
par.mlt = mlt;
par.mln = mln;
par.mln = mln;
par.mdx = mdx;
par.mdy = mdy;
par.mdz = mdz;

par.xx = xx;
par.yy = yy;
par.zz = zz;

par.nx = length(xx);
par.ny = length(yy);
par.nz = length(zz);
par.nmodel = par.nx*par.ny*par.nz;

par.min_x = min(xx);
par.max_x = max(xx);
par.min_y = min(yy);
par.max_y = max(yy);
par.min_z = min(zz);
par.max_z = max(zz); par.zmax = par.max_z; % not sure about adding this line, but weird to have two different vals

%% average velocities 
% calculated as mean over half-bins above and below par.zz

if par.PS == 1 || par.PS == 2 % just one velocity type
    par.vz = zeros(size(par.zz));
    zz = [0;par.zz;par.zz(end)];
    for iz = 1:length(par.vz)
        zbo = 0.5*(zz(iz+1)+zz(iz+2));
        zto = 0.5*(zz(iz)+zz(iz+1));
        par.vz(iz) = quad(@(z)vel_profile(par.PS,z),zto,zbo)/(zbo-zto);    
    end
    par.mvav   = interp1(par.zz,par.vz,par.mz);
elseif par.PS==3
    par.vz = zeros(length(par.zz),2);
    zz = [0;par.zz;par.zz(end)];
    for ip = 1:2
        for iz = 1:length(par.vz)
            zbo = 0.5*(zz(iz+1)+zz(iz+2));
            zto = 0.5*(zz(iz)+zz(iz+1));
            par.vz(iz,ip) = quad(@(z)vel_profile(ip,z),zto,zbo)/(zbo-zto);    
        end
        par.mvav(:,ip)   = interp1(par.zz,par.vz(:,ip),par.mz);
    end
    
end

%% average q 
% calculated as mean over half-bins above and below par.zz
Qmod = ql6;

Qmu = linterp(Qmod.depth,Qmod.qu,par.mz);
Qka = linterp(Qmod.depth,Qmod.qk,par.mz);
Qka(isinf(Qka)) = 9999; % reset to non-inf;

L = (4/3)*(1./1.732).^2; % this gives a Qp/Qs of 2.25

if par.PS == 1
    par.mQav = 1./(L./Qmu + (1-L)./Qka);
elseif par.PS == 2
    par.mQav = Qmu;
elseif par.PS == 3
    if par.Rdvpdvs == 0 % solving for both Qs and Qp
        par.mQav(:,1) = 1./(L./Qmu + (1-L)./Qka);
        par.mQav(:,2) = Qmu;
    else % scaling Qp, only solve for Qs
        par.mQav = Qmu;
    end
end

% [~,par] = make_start_model(par); % 

end
% 
% function V = vel_profile_integrand(z)
%     V = vel_profile(par.PS,z);
% end
