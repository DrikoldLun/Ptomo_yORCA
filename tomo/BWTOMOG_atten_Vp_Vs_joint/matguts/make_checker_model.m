function [ checker_model ] = make_checker_model(par,dval)


mf = zeros(par.ny,par.nx,par.nz);
k = -0.5;
gap = 100*sqrt(k^2+1); % between different bands
width = 50*sqrt(k^2+1); % within 1 band
b1 = 50; % y=x cutoff
b2 = -700; % y=x cutoff
b3 = 770; % y=kx cutoff
b4 = 200; % y=kx cutoff 300
mfxy = zeros(par.ny,par.nx);
for iy = 1:par.ny
y = par.yy(iy);
for ix = 1:par.nx
x = par.xx(ix);
op = (mod(floor((y-k*x)/gap),2)-0.5)*2; %anomaly polarity
if (mod(y-k*x,gap) <= width) && (-1/k*x+b2-y<=0) && (-1/k*x+b1-y>=0) && (b3+k*x-y>=0) && (b4+k*x-y<=0)
mfxy(iy,ix) = op;
end
end
end


%% vertical spacing
%mfz = [ 0  1  1  1  0  0 -1 -1  -1  0 0 1 1 1 ];
%mfz = [1 1 0 1 1 0];
mfz = [0 0 -1 -1 0 0 1 1 0];
%mfz = [ 0  0  1  1  1  0 0 -1  -1  -1 0 0 0 0 ];
%mfz = [ 0  0  1  1  1  0 0 0  0  0 0 0 0 0 ];
mfz = repmat(mfz,1,10);
%mfz = cat(mfz,flip(mfz));
mfz = mfz(1:par.nz);

%% now rep over layers
for iz = 1:par.nz
    mf(:,:,iz) = mfxy*mfz(iz);
end

mval = (0.01*dval)*mf(:); % 0.01 to convert from percent


%add a little random noise
% mval = mval + random('norm',0,0.00001,size(mval));

checker_model.mval = mval;
   
