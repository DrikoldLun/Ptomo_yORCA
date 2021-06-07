function [ H_smooth ] = make_smoothmat_dep( par,dep,direction )
% smooth using first derivative!!:

fprintf('>  creating smoothing matrix...\n')

zvh = par.smth_zvh; %0.3 ratio of vertical to horizontal smoothing

%% smooth nodes horizontally
%% xdirection
nnx = (par.nx-1)*par.ny*par.nz;
%rxi = zeros(2*nnx,1);
%rxj = zeros(2*nnx,1);
%rx  = zeros(2*nnx,1);
k = 0;
for ix = 1:par.nx-1
    for iy = 1:par.ny
        for iz = 1:par.nz
	    if (strcmp(direction,'down') && (par.zz(iz) > dep))||(strcmp(direction,'up') && (par.zz(iz) < dep))||(strcmp(direction,'both') && (par.zz(iz) < dep(1) || par.zz(iz) > dep(2)))
                continue;
	    else
                k = k+1;
                indx = find(par.mx==par.xx(ix) & par.my==par.yy(iy) & par.mz==par.zz(iz));
                indxp1 = find(par.mx==par.xx(ix+1) & par.my==par.yy(iy) & par.mz==par.zz(iz));
                rxi(2*k + [-1 0]) = k; % put two into k-th row
                rxj(2*k + [-1 0]) = [indx indxp1]; % put in column related to ix and ix+1
                rx( 2*k + [-1 0]) = [1 -1]*par.dh/abs(par.xx(ix+1)-par.xx(ix)); % weight by 1/distance between pt ix and ix+1
	    end
        end
    end
end

%% ydirection
nny = par.nx*(par.ny-1)*par.nz;
%ryi = zeros(2*nny,1);
%ryj = zeros(2*nny,1);
%ry  = zeros(2*nny,1);
k = 0;
for ix = 1:par.nx
    for iy = 1:par.ny-1
        for iz = 1:par.nz
	    if (strcmp(direction,'down') && (par.zz(iz) > dep))||(strcmp(direction,'up') && (par.zz(iz) < dep))||(strcmp(direction,'both') && (par.zz(iz) < dep(1) || par.zz(iz) > dep(2)))
                continue;
	    else
                k = k+1;
                indy = find(par.mx==par.xx(ix) & par.my==par.yy(iy) & par.mz==par.zz(iz));
                indyp1 = find(par.mx==par.xx(ix) & par.my==par.yy(iy+1) & par.mz==par.zz(iz));
                ryi(2*k + [-1 0]) = k; % put two into k-th row
                ryj(2*k + [-1 0]) = [indy indyp1]; % put in column related to iy and iy+1
                ry( 2*k + [-1 0]) = [1 -1]*par.dh/abs(par.yy(iy+1)-par.yy(iy)); % weight by 1/distance between pt iy and iy+1
	    end
        end
    end
end

%% smooth nodes vertically
%% zdirection
nnz = par.nx*par.ny*(par.nz-1);
%rzi = zeros(2*nnz,1);
%rzj = zeros(2*nnz,1);
%rz  = zeros(2*nnz,1);
k = 0;
for ix = 1:par.nx
    for iy = 1:par.ny
        for iz = 1:par.nz-1
	    if (strcmp(direction,'down') && (par.zz(iz+1) > dep))||(strcmp(direction,'up') && (par.zz(iz) < dep))||(strcmp(direction,'both') && (par.zz(iz) < dep(1) || par.zz(iz+1) > dep(2)))
		%disp(sprintf('%d,%d',par.zz(iz+1),dep))
		continue;
		%rz( 2*k + [-1 0]) = [0 0];
	    else
		k = k+1;
		indz = find(par.mx==par.xx(ix) & par.my==par.yy(iy) & par.mz==par.zz(iz));
		indzp1 = find(par.mx==par.xx(ix) & par.my==par.yy(iy) & par.mz==par.zz(iz+1));
		rzi(2*k + [-1 0]) = k; % put two into k-th row
            	rzj(2*k + [-1 0]) = [indz indzp1]; % put in column related to iz and iz+1
                rz( 2*k + [-1 0]) = [1 -1]*par.dh/abs(par.zz(iz+1)-par.zz(iz)); % weight by 1/distance between pt iz and iz+1
	    end
        end
    end
end


%% Compile
Rx = sparse(rxi,rxj,rx,length(unique(rxi)),par.nmodel);
Ry = sparse(ryi,ryj,ry,length(unique(ryi)),par.nmodel);

if k > 0
    Rz = sparse(rzi,rzj,rz,length(unique(rzi)),par.nmodel);
    H_smooth = [ par.dampx*Rx ;       % smooths in x-dir
             	 Ry ;       % smooths in y-dir
                 Rz*zvh ];  % smooths in z-dir
else
    H_smooth = [ par.dampx*Rx ;       % smooths in x-dir
                 Ry];       % smooths in y-dir
end
end
