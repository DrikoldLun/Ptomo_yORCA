% function [ vfield, rmserror ] = smooth_vctr_fld_invert( vx, vy, vth, vr, thopt, xnodes, ynodes, plotopt )
% % [ vfield ] = smooth_vctr_fld_invert( vx, vy, vth, vr, thopt, xnodes, ynodes, plotopt )
% % 
% % For a number of input vectors within a region, solve using a weighted,
% % damped version of Newton's method to arrive at an answer for the mixed
% % determined inverse problem. 
% %
% % Use formulation d=g(m) where the bottom part of g(m) is the prior
% % constraint that the change between adjacent model parameters is small,
% % weighted in such a way as to give epsilon weighting to
% % horizontal/vertical differences and epsilon/sqrt(2) to diagonals
% %
% %
% % INPUTS
% % vx = x-coordinates of vector centres
% % vy = y-coordinates of vector centres
% % vth = angles (theta) of vectors from "North" (the y-axis), in degrees
% % vr = magnitudes of vectors
% % thopt = option for vector directions 'direct' or 'axial'
% % xnodes = x-nodes for average grid
% % ynodes = y-nodes for average grid
% % plotopt = 0 (default) will not plot, 1 will make plot
% % 
% % OUTPUT
% % vfield = results matrix of average vectors with three layers
% %           - L1 is coordinates: x axis is real, y axis is imag
% %           - L2 has average vector directions
% %           - L3 has average vector lengths
% % rmserror = circular standard deviations of vectors in each averaging box

%% While testing
nobs=10;
vx = random('unif',0,2,nobs,1);
vy = random('unif',0,2,nobs,1);
vth = random('norm',90,2,nobs,1);
vr = random('unif',0,1,nobs,1);
xnodes = 0:1:2;
ynodes = 0:1:2;
thopt = 'direct';
plotopt = 0;
%%

if nargin < 8 
    plotopt = 0;
end


if length(vx)==length(vy) && length(vx)==length(vth) && length(vx)==length(vr)
    if max(vx)<=max(xnodes) && min(vx)>=min(xnodes)...
            && max(vy)<=max(ynodes) && min(vy)>=min(ynodes)
        fprintf('Good to go\n')
    else
        fprintf('N.B. some vectors are out of node bounds and will be ignored\n')
    end
else
    error('Some inputs are wrong size')
end

if strcmp(thopt,'axial')==1
    vth = mod(2.*vth,360);
end


n = length(vx); N = 2*n;
nxn=length(xnodes); X = nxn-1;
nyn=length(ynodes); Y = nyn-1;
rmax=max(vr);

M = 2*X*Y;

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
% X = 3;
% Y = 3;
% n = 3;

vfield = zeros(Y,X,2);
vfield(:,:,1) = magic(X); vfield(:,:,2) = random('unid',10,Y,X); % need some starting values...
epsilon = 1;
M = 2*X*Y;
NN = 2*n + 2*(4*X*Y - 3*X - 3*Y + 2);
%% TESTING

%% set up data vector Cs then Ss - rest is differences, so zeros
d = zeros(NN,1);
for in = 1:n
    d(in) = vr(in)*cosd(vth(in)); % first half is cos
    d(in + NN/2) = vr(in)*sind(vth(in)); % second half is sin
end

mest = zeros(M,1);
for ix = 1:X
    for iy = 1:Y
        mest(iy + (ix-1)*Y) = vfield(Y+1-iy,ix,1);
        mest(M/2 + iy + (ix-1)*Y) = vfield(Y+1-iy,ix,2);
    end
end

E = 0;
dE = 1000;
%% gm
% gm is thus a vector, NN elements long. the first NN/2 elements are the
% predicted phis, followed by the difference in phis, and then the same for
% dts
gm = zeros(NN,1);
%% differences in g(m)
% diffs up
for ix = 1:X
    for iy = 1:Y-1
        kk = (ix-1)*Y + iy; % position in m
        KK = n + (ix-1)*(Y-1) + iy;        % position in gm
        gm(KK) = epsilon*(mest(kk) - mest(kk + 1)); % phis
        gm(NN/2 + KK) = epsilon*(mest(M/2 + kk) - mest(M/2 + kk + 1)); % dts
    end
end
% diffs right
for ix = 1:X-1
    for iy = 1:Y
        kk = (ix-1)*Y + iy; % position in m
        KK = n + X*(Y-1) + (ix-1)*Y + iy;        % position in gm
        gm(KK) = epsilon*(mest(kk) - mest(kk + Y)); % phis
        gm(NN/2 + KK) = epsilon*(mest(M/2 + kk) - mest(M/2 + kk + Y)); % dts    
    end
end
% diffs diag up&right
for ix = 1:X-1
    for iy = 1:Y-1
        kk = (ix-1)*Y + iy; % position in m
        KK = n + X*(Y-1) + (X-1)*Y + (ix-1)*(Y-1) + iy;        % position in gm
        gm(KK) = (epsilon/sqrt(2))*(mest(kk) - mest(kk + Y + 1));
        gm(NN/2 + KK) = (epsilon/sqrt(2))*(mest(M/2 + kk) - mest(M/2 + kk + Y +1));
    end
end
% diffs diag up&left
for ix = 2:X
    for iy = 1:Y-1
        kk = (ix-1)*Y + iy; % position in m
        KK = n + X*(Y-1) + (X-1)*Y + (X-1)*(Y-1) + (ix-2)*(Y-1) + iy;        % position in gm
        gm(KK) = (epsilon/sqrt(2))*(mest(kk) - mest(kk - Y + 1));
        gm(NN/2 + KK) = (epsilon/sqrt(2))*(mest(M/2 + kk) - mest(M/2 + kk - Y +1));
    end
end
% dest
for in = 1:n % loop over data
gm(in) = mest(boxij(in) + M/2)* cosd(mest(boxij(in)));
gm(in + NN/2) = mest(boxij(in) + M/2)* sind(mest(boxij(in)));
end

vfield
[gm(1:NN/2) gm(NN/2 + 1:end)]

%% Gm
% Gm is a matrix, NN x M. the first NN/2 rows are the
% predicted phis, followed by the difference in phis, and then the same for
% dts. Each column is the element of gm differentiated by m_j
Gm = spalloc(NN,M,2*NN);
Gm = zeros(NN,M);
%% differences in g(m)
% diffs up
for ix = 1:X
    for iy = 1:Y-1
        kk = (ix-1)*Y + iy; % position in m
        KK = n + (ix-1)*(Y-1) + iy;        % position in Gm
        Gm(KK,kk) = epsilon; % phis
        Gm(KK,kk + 1) = -epsilon; % phis
        Gm(NN/2 + KK,M/2 + kk) = epsilon; % dts
        Gm(NN/2 + KK,M/2 + kk + 1) = - epsilon; % dts
    end
end
% diffs right
for ix = 1:X-1
    for iy = 1:Y
        kk = (ix-1)*Y + iy; % position in m
        KK = n + X*(Y-1) + (ix-1)*Y + iy;        % position in gm
        Gm(KK,kk) = epsilon; % phis
        Gm(KK,kk + Y) = -epsilon; % phis
        Gm(NN/2 + KK,M/2 + kk) = epsilon; % dts    
        Gm(NN/2 + KK,M/2 + kk + Y) = -epsilon; % dts    
    end
end
% diffs diag up&right
for ix = 1:X-1
    for iy = 1:Y-1
        kk = (ix-1)*Y + iy; % position in m
        KK = n + X*(Y-1) + (X-1)*Y + (ix-1)*(Y-1) + iy;        % position in gm
        Gm(KK,kk) = epsilon/sqrt(2);
        Gm(KK,kk + Y + 1) = -epsilon/sqrt(2);
        Gm(NN/2 + KK,M/2 + kk) = epsilon/sqrt(2);
        Gm(NN/2 + KK,M/2 + kk + Y +1) = -epsilon/sqrt(2);
    end
end
% diffs diag up&left
for ix = 2:X
    for iy = 1:Y-1
        kk = (ix-1)*Y + iy; % position in m
        KK = n + X*(Y-1) + (X-1)*Y + (X-1)*(Y-1) + (ix-2)*(Y-1) + iy;        % position in gm
        Gm(KK,kk) = epsilon/sqrt(2);
        Gm(KK,kk - Y + 1) = -epsilon/sqrt(2);
        Gm(NN/2 + KK,M/2 + kk) = epsilon/sqrt(2);
        Gm(NN/2 + KK,M/2 + kk - Y +1) = -epsilon/sqrt(2);
    end
end
% dest
for in = 1:n % loop over data
Gm(in,boxij(in) + M/2) = cosd(mest(boxij(in)));
Gm(in,boxij(in)) = mest(boxij(in) + M/2)* -sind(mest(boxij(in)));
Gm(in + NN/2,boxij(in) + M/2) = sind(mest(boxij(in)));
Gm(in + NN/2,boxij(in)) = mest(boxij(in) + M/2)* cosd(mest(boxij(in)));
end

% %% TESTING
% 
% for ix = 1:nxn-1
%     for iy = 1:nyn-1
%         vfield(iy,ix,1) = 1i*mean([ynodes(iy),ynodes(iy+1)]) + mean([xnodes(ix),xnodes(ix+1)]);
%         
% %         z = intersect(find(vx<xnodes(ix+1) &  vx>xnodes(ix)),find(vy<xnodes(iy+1) &  vy>xnodes(iy)));
%         z = intersect(find(sqrt((vx-xnodes(ix+1)).^2 + (vx-xnodes(ix)).^2) < mean(diff(xnodes))),...
%                       find(sqrt((vy-ynodes(iy+1)).^2 + (vy-ynodes(iy)).^2) < mean(diff(ynodes))));
% 
%         if isempty(z)==1; continue; end
%         zr=vr(z);
%         zth=vth(z);
% %         if strcmp(thopt,'axial')==1
% %             zth = mod(vth(z)*2,360);
% %         elseif strcmp(thopt,'direct')==1
% %             zth = vth(z);
% %         end
%         
%         [vfield(iy,ix,3),vfield(iy,ix,2),rmserror(iy,ix)] = mrl(d2r(zth),zr,thopt);
%         
%         if strcmp(thopt,'axial')==1
%             vfield(iy,ix,2) = r2d(vfield(iy,ix,2));
%         elseif strcmp(thopt,'direct')==1
%             vfield(iy,ix,2) = r2d(vfield(iy,ix,2));
%         end
%         
%     end
% end

if plotopt==1
figure(1); clf
for ii=1:n
    x1 = vr(ii)*sind(vth(ii)) + vx(ii); y1 = vr(ii)*cosd(vth(ii)) + vy(ii);
    x2 =-vr(ii)*sind(vth(ii)) + vx(ii); y2 =-vr(ii)*cosd(vth(ii)) + vy(ii);
    hold on;	line([x1 x2],[y1 y2],'linewidth',1,'color','r');	hold off
end
set(gca,'XLim',[xnodes(1)-rmax xnodes(end)+rmax],'YLim',[ynodes(1)-rmax ynodes(end)+rmax],'XTick',xnodes,'YTick',ynodes)
end

if plotopt==1
for ix = 1:nxn-1
    for iy = 1:nyn-1
        if vfield(iy,ix,3)==0, continue; end
    x1 =  vfield(iy,ix,3)*sind(vfield(iy,ix,2)) + real(vfield(iy,ix,1)); y1 =  vfield(iy,ix,3)*cosd(vfield(iy,ix,2)) + imag(vfield(iy,ix,1));
    x2 = -vfield(iy,ix,3)*sind(vfield(iy,ix,2)) + real(vfield(iy,ix,1)); y2 = -vfield(iy,ix,3)*cosd(vfield(iy,ix,2)) + imag(vfield(iy,ix,1));
    hold on;	line([x1 x2],[y1 y2],'linewidth',2,'color','b');	hold off
    end
end
end




% end

