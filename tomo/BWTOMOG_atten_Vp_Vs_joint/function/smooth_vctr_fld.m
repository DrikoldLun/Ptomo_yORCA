function [ vfield, rmserror ] = smooth_vctr_fld( vx, vy, vth, vr, thopt, xnodes, ynodes, plotopt )
% [ vfield ] = smooth_vctr_fld( vx, vy, vth, vr, thopt, xnodes, ynodes, plotopt )
% 
% For a number of input vectors within a region, divide them up into boxes
% and do a sort of running average vector over each unit of the subdivided
% space
%
% INPUTS
% vx = x-coordinates of vector centres
% vy = y-coordinates of vector centres
% vth = angles (theta) from "North" (the y-axis) in degrees of vectors
% vr = magnitudes of vectors
% thopt = option for vector directions 'direct' or 'axial'
% xnodes = x-nodes for average grid
% ynodes = y-nodes for average grid
% plotopt = 0 (default) will not plot, 1 will make plot
% 
% OUTPUT
% vfield = results matrix of average vectors with three layers
%           - L1 is coordinates: x axis is real, y axis is imag
%           - L2 has average vector directions
%           - L3 has average vector lengths
% rmserror = circular standard deviations of vectors in each averaging box


% %% While testing
% vx = random('unif',0,10,800,1);
% vy = random('unif',0,10,800,1);
% vth = random('norm',90,10,800,1);
% vr = random('unif',0,1,800,1);
% xnodes = 0:1:10;
% ynodes = 0:1:10;
% thopt = 'direct';
% plotopt = 1;
% %%

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

n = length(vx);
nxn=length(xnodes);
nyn=length(ynodes);
rmax=max(vr);

if plotopt==1
figure(1); clf
for ii=1:n
    x1 = vr(ii)*sind(vth(ii)) + vx(ii); y1 = vr(ii)*cosd(vth(ii)) + vy(ii);
    x2 =-vr(ii)*sind(vth(ii)) + vx(ii); y2 =-vr(ii)*cosd(vth(ii)) + vy(ii);
    hold on;	line([x1 x2],[y1 y2],'linewidth',1,'color','r');	hold off
end
set(gca,'XLim',[xnodes(1)-rmax xnodes(end)+rmax],'YLim',[ynodes(1)-rmax ynodes(end)+rmax],'XTick',xnodes,'YTick',ynodes)
end

% results structure - the first layer is the coordinates (complex) 
% the second has the averaged vector angles, third has the vector mrls
vfield = zeros(nyn-1,nxn-1,3); 
rmserror = zeros(nyn-1,nxn-1); 
for ix = 1:nxn-1
    for iy = 1:nyn-1
        vfield(iy,ix,1) = 1i*mean([ynodes(iy),ynodes(iy+1)]) + mean([xnodes(ix),xnodes(ix+1)]);
        
%         z = intersect(find(vx<xnodes(ix+1) &  vx>xnodes(ix)),find(vy<xnodes(iy+1) &  vy>xnodes(iy)));
        z = intersect(find(sqrt((vx-xnodes(ix+1)).^2 + (vx-xnodes(ix)).^2) < mean(diff(xnodes))),...
                      find(sqrt((vy-ynodes(iy+1)).^2 + (vy-ynodes(iy)).^2) < mean(diff(ynodes))));

        if isempty(z)==1; continue; end
        zr=vr(z);
        zth=vth(z);
%         if strcmp(thopt,'axial')==1
%             zth = mod(vth(z)*2,360);
%         elseif strcmp(thopt,'direct')==1
%             zth = vth(z);
%         end
        
        [vfield(iy,ix,3),vfield(iy,ix,2),rmserror(iy,ix)] = mrl(d2r(zth),zr,thopt);
        
        if strcmp(thopt,'axial')==1
            vfield(iy,ix,2) = r2d(vfield(iy,ix,2));
        elseif strcmp(thopt,'direct')==1
            vfield(iy,ix,2) = r2d(vfield(iy,ix,2));
        end
        
    end
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




end

