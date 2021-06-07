function [ H_damp,h_damp ] = make_dampmat_dep( data,par,dep,direction )
% makes damping matrix 
% push layer average to zero
% push static terms to zero
% 
% model parameter vector is [nmodel,nevts,nstas]

fprintf('>  Creating damping matrix... \n');

damp_lav0 = 10; % damping layer average to zero
damp_age = 100; % smoothing along age contours
damp_zlt40 = 10; % damp shallow most nodes more

% account for doubling of station terms if P and S
if par.PS == 3, sdim = 2; else, sdim = 1; end
if par.PS == 3 && par.Rdvpdvs==0, mdim = 2; else mdim = 1; end

nstas = data.stn.nstas*sdim;   
nevts = data.evt.nevts;
nmod = par.nmodel*mdim;                          
nmat = nmod+nevts+nstas; % horizontal dimension of smoothing matrix

zz = par.zz;

%% OVERALL damping (i.e. x,y=0)
% make weights for diagonal - damping all to zeros
mwt = ones(nmod,1);
zuse = repmat(par.mz,mdim,1); % double up length vector if needed
for iz = 1:size(par.depsmth,1)
    mwt((zuse >= par.depsmth(iz,1)) & (zuse <= par.depsmth(iz,2))) = par.depsmth(iz,3);
end

%{
if strcmp(direction,'down')
    mwt(par.mz>dep) = 0;
elseif strcmp(direction,'up')
    mwt(par.mz<dep) = 0;
end
%}
if strcmp(direction,'down')
    ind = find(par.mz<=dep);
elseif strcmp(direction,'up')
    ind = find(par.mz>=dep);
elseif strcmp(direction,'both')
    ind = find(par.mz>=dep(1)&par.mz<=dep(2));
end

% mwt(par.mz<60) = 2; % extra damp for shallow nodes
% mwt(par.mz>350) = 2; % extra damp for deep nodes

% turn into sparse diagonal matrix
for i = 1:length(ind)
    si(i) = i;
    sj(i) = ind(i);
    s(i) = par.damp*mwt(ind(i));
end
%
%si = 1:nmod;
%sj = 1:nmod;
%s  = par.damp*mwt;
%
% fprintf('     NB - NO DIAGONAL DAMPING\n')

H_damp = sparse(si,sj,s,length(ind),nmat);
h_damp = zeros(length(ind),1);

% %% Shallow node damping
% shindx = find(par.mz<40);
% N = length(shindx);
% mwt = ones(N,1);
% 
% shi = shindx;
% shj = 1:N;
% sh  = damp_zlt40*mwt;
% 
% shD = sparse(shi,shj,sh,N,nmat);
% 
% H_damp = [H_damp;   shD                      ];
% h_damp = [h_damp;   zeros(N,1) ];

%% add event damping 
ewt = ones(nevts,1);

eadi = [1:nevts]';
eadj = [1:nevts]'+nmod;
ead  = par.damp_evt*ewt;
ED = sparse(eadi,eadj,ead,nevts,nmat);

H_damp = [H_damp;   ED                      ];
h_damp = [h_damp;   par.start_model.estatic ];


%% add station damping 
swt = ones(nstas,1);
sadi = [(1:nstas)]';
sadj = [(1:nstas)]'+nmod+nevts;
sad  = par.damp_stn*swt;
SD = sparse(sadi,sadj,sad,nstas,nmat);

H_damp = [H_damp;   SD                      ];
h_damp = [h_damp;   par.start_model.sstatic ];

%% add constraint that average in each layer is 0
bi = zeros(nmod,1);
bj = zeros(nmod,1);
b  = zeros(nmod,1);
nlay = par.nx*par.ny;
n = 1;
for ip = 1:mdim
for iz = 1:par.nz 
    bi(n:n+nlay-1) = (iz + (ip-1)*par.nz) * ones(nlay,1);
    bj(n:n+nlay-1) = find(par.mz==zz(iz)) + (ip-1)*par.nmodel;
    b(n:n+nlay-1)  = ones(nlay,1);
    n = n+nlay;
end
end
B = sparse(bi,bj,b,par.nz*mdim,nmat);

H_damp = [H_damp;   damp_lav0*B     ];
h_damp = [h_damp;   zeros(par.nz*mdim,1) ];
    
%% Smooth by ocean floor age
% smooth in age bins of 2 Ma. Only do for ages of < 12 Ma
if par.age_smooth
disp('     smoothing along isochrons ');
    
ages = 0:2:12;
Nag = length(ages)-1;
ci = nan(2*nmod,1); % just to be safe, max #�of elements is 2*nmod
cj = nan(2*nmod,1);
c  = nan(2*nmod,1);

n = 1; % number of elements in sparse index vectors
iazi = 0; % number of age-depth perms
for iz = 1:par.nz
for ia = 1:Nag
    % find points in this depth slice within the age bounds
    inds = find(par.mage <= ages(ia+1) & par.mage >= ages(ia) & par.mz == zz(iz)); 
    ninds = length(inds);
    inds = [inds(:);inds(1)]; % just tack on first element at the end to complete loop

    %     Now string these all together by looping through them minimising
    %     first difference between each and next
    for ii = 1:ninds
        iazi = iazi+1;
        
        ci(n+[0:1]) = iazi*[1 1];
        cj(n+[0:1]) = [inds(ii) inds(ii+1)];
        c(n+[0:1])  = [-1 1] ;
        
        n = n+2;
    end
end
end
ci(isnan(ci))=[];
cj(isnan(cj))=[];
c(isnan(c))=[];

C = sparse(ci,cj,c,iazi,nmat);

H_damp = [H_damp;   damp_age*C      ];
h_damp = [h_damp;   zeros(iazi,1)   ];

end


end

    

