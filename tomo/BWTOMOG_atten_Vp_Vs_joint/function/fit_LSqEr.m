function [ m,b,mstd,bstd,m_upper,m_lower ] = fit_LSqEr( X,Y,b0opt,normopt,mprior,bprior,xwt,ywt )
%[ m,b,mstd,bstd,m_upper,m_lower ] = fit_LSqEr( X,Y,b0opt=1,normopt=2,mprior,bprior,xwt,ywt )
%   function to fit a straight line to X,Y data where there is error in
%   both data by grid searching to find the line that gives the minimum
%   squared misfit in both dimensions with a line equation y = m*x + b
%
% Expand in future to do F-test on grid search to estimate uncertainties to
% the fit, and to add options for L1 norm instead of L2 norm

if nargin<3 || isempty(b0opt)
    b0opt = 1; %option to pin b to zero
end
if nargin<4 || isempty(normopt)
    normopt = 2;
end
if nargin < 5 || isempty(mprior)
    % find something about the bounds of the data
    dx = max(X)-min(X);
    dy = max(Y)-min(Y);
    m0 = dy/dx;
    mprior = -1000:1000;
end
if nargin < 6 || isempty(bprior)
    if b0opt
        bprior = 0;
    else
        bprior = -1000:1000;
    end
end
if nargin < 7
    xwt = ones(size(X));
end
if nargin < 8
    ywt = ones(size(X));
end

M = length(mprior);
B = length(bprior);

misfit = zeros(M,B);

for im = 1:length(mprior)
    for ib = 1:length(bprior)
        dy = Y - ( mprior(im)*X + bprior(ib) );
        dx = X - ( (Y-bprior(ib))./mprior(im) );
        dy = ywt.*dy;
        dx = xwt.*dx;
        misfit(im,ib) = norm(dy,normopt) + norm(dx,normopt);
    end
end

[misfit_best,ib_best,im_best] = mingrid(misfit);

m = mprior(im_best);
b = bprior(ib_best);

%% F tests to get confidence bounds for m
nar = length(X);
Fval = zeros(M,1);
for ii = 1:M
    Fval(ii) = (misfit(im_best,ib_best)\misfit(ii,ib_best));
end
Fcrit = finv(0.67,nar-1,nar-1);
m_upper = mprior(find(Fval<Fcrit,1,'last'));
m_lower = mprior(find(Fval<Fcrit,1,'first'));

%% F tests to get confidence bounds for b
nar = length(X);
Fval = zeros(B,1);
for ii = 1:B
    Fval(ii) = (misfit(im_best,ib_best)\misfit(im_best,ii));
end
Fcrit = finv(0.67,nar-1,nar-1);
b_upper = bprior(find(Fval<Fcrit,1,'last'));
b_lower = bprior(find(Fval<Fcrit,1,'first'));


mstd = (m_upper-m_lower)/2;
bstd = (b_upper-b_lower)/2;

end

