function [ vfield, rmserror ] = smooth_vctr_fld_invert( vx, vy, Emap, xnodes, ynodes, plotopt )
% [ vfield ] = smooth_vctr_fld_invert( vx, vy, Emap, xnodes, ynodes, plotopt )
% 
% For a number of input vectors within a region, solve using a weighted,
% damped version of Newton's method to arrive at an answer for the mixed
% determined inverse problem. 
%
% Use formulation d=g(m) where the bottom part of g(m) is the prior
% constraint that the change between adjacent model parameters is small,
% weighted in such a way as to give epsilon weighting to
% horizontal/vertical differences and epsilon/sqrt(2) to diagonals
%
% NB the function seeks to minimise the values on the Emap, assuming they
% are errors - for RC method, flip to make max correlation the minimum
% value
% NB assume max dt is 4 seconds
%
%
% INPUTS
% vx = x-coordinates of vector centres
% vy = y-coordinates of vector centres
% Emap = maps of error with phi and dt - a matrix [M_phi,M_dt,N_obs]
% xnodes = x-nodes for average grid
% ynodes = y-nodes for average grid
% plotopt = 0 (default) will not plot, 1 will make plot
% 
% OUTPUT
% vfield = results matrix of average vectors with three layers
%           - L1 is coordinates: x axis is real, y axis is imag
%           - L2 has average vector directions
%           - L3 has average vector lengths
% rmserror = rms error summed from all inputs

%% While testing
% nobs=50;
% xlims = [0,3];
% ylims = [0,3];
% vx = random('unif',xlims(1),xlims(2),nobs,1);
% vy = random('unif',ylims(1),ylims(2),nobs,1);
% x_th = -90:2:90;
% x_r = 0:0.04:4;
% 
% xnodes = xlims(1):1:xlims(2);
% ynodes = ylims(1):1:ylims(2);
% 
% Emap = zeros(length(x_th),length(x_r),nobs);
% dt_tru = zeros(nobs,1);
% phi_tru = zeros(nobs,1);
% for in = 1:nobs
% vth = random('norm',vx(in)*45,5); phi_tru(in) = vth;
% vr = random('norm',2,0.3);        dt_tru(in) = vr;
% [Emap(:,:,in),EX,EY] = gauss2D(vr,0.8,vth,15,x_r,x_th);
% % surf(x_r,x_th,Emap(:,:,in))
% % pause(5)
% end
% Emap = 1-Emap;
% vx = ilons;
% vy = ilats;
% Emap = ESCmaps;
% xnodes = [lonlim(1):gridspacing:lonlim(2)];
% ynodes = [latlim(1):gridspacing:latlim(2)];
% 
% 
% 
% nargin=5;
% plotopt = 1;
%% True start


if nargin < 5 
    plotopt = 0;
end

% check lengths of vectors match and points are inside grid
if length(vx)==length(vy) && length(vx)==size(Emap,3) 
    goodin = find(vx>=min(xnodes)&vx<=max(xnodes)&vy>=min(ynodes)&vy<=max(ynodes));
    if length(goodin)<length(vx)
        fprintf('N.B. some vectors are out of node bounds and will be ignored\n')
        vx = vx(goodin); 
        vy = vy(goodin);
        Emap = Emap(:,:,goodin);
    else
        fprintf('Good to go\n')
    end
else
    error('Some inputs are wrong size')
end

% shape x and y node vectors
xnodes = reshape(xnodes,1,length(xnodes)); % row vector
ynodes = reshape(ynodes,length(ynodes),1); % column vector

%values corresponding to each row and column of Emap
if isodd(size(Emap,1))==1
x_phi = flipud([-90:180/(size(Emap,1)-1):90]');
else
x_phi = flipud([-90+180/(size(Emap,1)):180/(size(Emap,1)):90]');
end
x_dt = [0:4/(size(Emap,2)-1):4];


n = length(vx); N = 2*n;
nxn=length(xnodes); X = nxn-1;
nyn=length(ynodes); Y = nyn-1;
rmax=2; % notional max length of vector - just to make plots nice

M = 2*X*Y; % number of model parameters

%% Use vx and vy to assoc. each obs with a given box
boxij=zeros(n,1);
for in = 1:n
    ix = find(xnodes < vx(in),1,'last');
    iy = find(ynodes < vy(in),1,'last');
    boxij(in) = iy + (ix-1)*Y;
end

% structure model parameter column vector such that it has successive
% columns from the spatial matrix (i.e. [c1;c2;c3;...]) going from the
% BOTTOM UPWARDS and has all the phis and then all the dTs
% Thus m = [phi11;phi21;phi31;...phiY1;phi12;phi22;...phiYX;dt11;dt21;dt31;...dtY1;dt12;dt22;...dtYX;
%
% d is similar but has all the C's first, then all S's, but is in order of
% input vectors - random w.r.t. space.
%   C = dT*cos(phi)
%   S = dT*sin(phi)

% damping for slow change in model parameters: m1-m2 of adjacent cells
% assumed to be zero. weight by epsilon for truly adjacent boxes, and with
% epsilon/sqrt(2) for diagonals. Only do each difference once - do this by
% taking differences going up and right. 
% going up,         X*(Y-1) diffs
% going right,      (X-1)*Y diffs
% going up&right   (X-1)*(Y-1) diffs
% going up&left    (X-1)*(Y-1) diffs
% = XY - X + XY - Y + 2XY - 2X - 2Y + 2 
% total of 4*X*Y - 3*X - 3*Y + 2
ndamping = (4*X*Y - 3*X - 3*Y + 2); % and NB this is for both phi and dt


% total length of data vector is thus:
NN = N + 2*ndamping;
dobs = zeros(NN,1);
% total length of model vector is thus:
mest = zeros(M,1);
% where mest(ix,iy) = mest(iy + (ix-1)*Y)
% i.e. mest(iy + (ix-1)*Y) = vfield(Y+1-iy,ix)
% Need some a priori values for the model parameters
%% !!!!!!!! see line above

% results structure - the first layer is the coordinates (complex) 
% the second has the averaged vector angles, third has the vector mrls
vfield = zeros(Y,X); 
rmserror = zeros(Y,X); 

%% TESTING

%% test parameters =========================
niter = 1000;
nhoming = 20; 
Emin = 1e11;
epsilon= 0.3e9;
phi2dt_scale = 4/180;
%% test parameters =========================

mestbest = zeros(Y,X,2);
for nh = 1:nhoming
fprintf('Homer is %u\n',nh)
for iter = 1:niter
vfield = zeros(Y,X,2);
if nh==1
vfield(:,:,1) = random('unif',-90,90,Y,X); % 1st layer is phi
vfield(:,:,2) = random('unif',0,4,Y,X);    % 2nd layer is dt
else
    % after the first set of iterations, take the starting estimates as the
    % previous best, and try some scatter about these estimates, with a
    % standard deviation proportional to the error in that box
    for ix = 1:X
        for iy =1:Y
            fact = Esmin(iy,ix)/Emin; if fact==0, fact=1; end
            vfield(iy,ix,1) = mestbest(iy,ix,1) + random('norm',0,90*fact); % 1st layer is phi
            vfield(iy,ix,2) = mestbest(iy,ix,2) + random('norm',0,2*fact);    % 2nd layer is dt
        end
    end
% vfield(:,:,1) = mestbest(:,:,1) + random('norm',0,20,Y,X); % 1st layer is phi
% vfield(:,:,2) = mestbest(:,:,2) + random('norm',0,1,Y,X);    % 2nd layer is dt    
% limit angles
vfield(:,:,1) = mod(vfield(:,:,1)+90,180)-90;
% junk = reshape(vfield(:,:,1),size(vfield(:,:,1),1)*size(vfield(:,:,1),2),1);
% junk(junk>90)=90;
% junk(junk<-90)=-90;
% vfield(:,:,1) = reshape(junk,size(vfield(:,:,1),1),size(vfield(:,:,1),2));
% limit times
% junk = reshape(vfield(:,:,2),size(vfield(:,:,2),1)*size(vfield(:,:,2),2),1);
% junk(junk>4)=4;
% junk(junk<0)=0;
% vfield(:,:,2) = reshape(junk,size(vfield(:,:,2),1),size(vfield(:,:,2),2));    
vfield(:,:,2) = mod(vfield(:,:,2),4);
end
    
% figure(2); hold on
% plot(vfield(:,:,1),vfield(:,:,2),'x')
% hold off
% ylim([-90,90]); xlim([0,4]);


M = 2*X*Y;
NN = 2*n + 2*(4*X*Y - 3*X - 3*Y + 2);
%% TESTING

mest = zeros(M,1);
for ix = 1:X
    for iy = 1:Y
        mest(iy + (ix-1)*Y) = vfield(Y+1-iy,ix,1);       % 1st M/2 is phis
        mest(M/2 + iy + (ix-1)*Y) = vfield(Y+1-iy,ix,2); % 2nd M/2 is dts
    end
end

E = 0;
dE = 1000;

%% differences in g(m)
% diffs up
du = zeros(2*X*(Y-1),1);
for ix = 1:X
    for iy = 1:Y-1
        kk = (ix-1)*Y + iy; % position in m
        KK = (ix-1)*(Y-1) + iy;        % position in du
        du(KK) = phi2dt_scale*epsilon*(mest(kk) - mest(kk + 1)); % phis
        du(X*(Y-1) + KK) = epsilon*(mest(M/2 + kk) - mest(M/2 + kk + 1)); % dts
    end
end
% diffs right
dr = zeros(2*(X-1)*Y,1);
for ix = 1:X-1
    for iy = 1:Y
        kk = (ix-1)*Y + iy; % position in m
        KK = (ix-1)*Y + iy;        % position in gm
        dr(KK) = phi2dt_scale*epsilon*(mest(kk) - mest(kk + Y)); % phis
        dr((X-1)*Y + KK) = epsilon*(mest(M/2 + kk) - mest(M/2 + kk + Y)); % dts    
    end
end
% diffs diag up&right
dur = zeros(2*(X-1)*(Y-1),1);
for ix = 1:X-1
    for iy = 1:Y-1
        kk = (ix-1)*Y + iy; % position in m
        KK = (ix-1)*(Y-1) + iy;        % position in gm
        dur(KK) = phi2dt_scale*(epsilon/sqrt(2))*(mest(kk) - mest(kk + Y + 1));
        dur((X-1)*(Y-1) + KK) = (epsilon/sqrt(2))*(mest(M/2 + kk) - mest(M/2 + kk + Y +1));
    end
end
% diffs diag up&left
dul = zeros(2*(X-1)*(Y-1),1);
for ix = 2:X
    for iy = 1:Y-1
        kk = (ix-1)*Y + iy; % position in m
        KK = (ix-2)*(Y-1) + iy;        % position in gm
        dul(KK) = phi2dt_scale*(epsilon/sqrt(2))*(mest(kk) - mest(kk - Y + 1));
        dul((X-1)*(Y-1) + KK) = (epsilon/sqrt(2))*(mest(M/2 + kk) - mest(M/2 + kk - Y +1));
    end
end

L = sum(abs([du;dr;dul;dur]))/(2*ndamping); % "length" of solution (actually epsilon*length/ndamped)

% sum of errors
Es = zeros(Y,X);
for in = 1:n % loop over data
    % put errors in grid matching results grid - i.e. error by box.
ix = 1 + floor((boxij(in)-1)/Y); % horiz position in grid
iy = 1 + mod(Y - boxij(in),Y); % vert position in grid (matrix convention, i.e. 1 is top row)
phiest = mest(boxij(in));
dtest = mest(boxij(in)+ M/2);
Es(iy,ix) = Emap(abs(x_phi--phiest)==min(abs(x_phi--phiest)),...
                abs(x_dt-dtest)==min(abs(x_dt-dtest)),in)/sum(boxij==boxij(in)); 
% error value is location on Emap corresponding to the estimated model params
% scale each error to 1 - divide by number of obs in that box when adding
end
E = sum(sum(Es)); % mean error


if E + L < Emin
    Esmin = Es;
    Emin = E+L
    mestbest = vfield;
end
end % loop on iters
end % loop on homing
% return

mestbest;
Esmin
Emin
vfield = zeros(Y,X,3);
vfield(:,:,2) = mestbest(:,:,1); % 2nd layer is phi
vfield(:,:,3) = mestbest(:,:,2); % 3rd layer is dt
[A,B] = meshgrid(midpts(xnodes),flipud(1i*midpts(ynodes)));
vfield(:,:,1) = A+B             % 1st layer is x-iy locationEmap

% if plotopt==1
% figure(1); clf
% scale = 0.25;
% for ii=1:n
%     x1 = scale*dt_tru(ii)*sind(phi_tru(ii)) + vx(ii); y1 = scale*dt_tru(ii)*cosd(phi_tru(ii)) + vy(ii);
%     x2 =-scale*dt_tru(ii)*sind(phi_tru(ii)) + vx(ii); y2 =-scale*dt_tru(ii)*cosd(phi_tru(ii)) + vy(ii);
%     hold on;	line([x1 x2],[y1 y2],'linewidth',1,'color','r');	hold off
% end
% set(gca,'XLim',[xnodes(1)-scale*rmax xnodes(end)+scale*rmax],'YLim',[ynodes(1)-scale*rmax ynodes(end)+scale*rmax],'XTick',xnodes,'YTick',ynodes)
% grid
% end


if plotopt==1
for ix = 1:nxn-1
    for iy = 1:nyn-1
        if vfield(iy,ix,3)==0, continue; end
    x1 =  scale*vfield(iy,ix,3)*sind(vfield(iy,ix,2)) + real(vfield(iy,ix,1)); y1 =  scale*vfield(iy,ix,3)*cosd(vfield(iy,ix,2)) + imag(vfield(iy,ix,1));
    x2 = -scale*vfield(iy,ix,3)*sind(vfield(iy,ix,2)) + real(vfield(iy,ix,1)); y2 = -scale*vfield(iy,ix,3)*cosd(vfield(iy,ix,2)) + imag(vfield(iy,ix,1));
    hold on;	line([x1 x2],[y1 y2],'linewidth',2,'color','b');	hold off
    end
end
end

% figure(3);
% pcolor(sum(Emap,3))
% hold on
% plot(find(abs(x_dt-mestbest(:,:,2))==min(abs(x_dt-mestbest(:,:,2)))),...
%     length(x_phi)-find(abs(x_phi-mestbest(:,:,1))==min(abs(x_phi-mestbest(:,:,1)))),'ow')
% hold off


% end

