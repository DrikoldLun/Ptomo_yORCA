function [ A, phi, rmserrorY, rmserrorX ] = fit_sinNtheta_L1( xx,yy,n,We,plotopt )
%function [ A, phi, rmserrorY, rmserrorX ] = fit_sinNtheta_L1( xx,yy,n,We,plotopt )
%OR 
%function [ A, phi, rmserrorY, rmserrorX ] = fit_sinNtheta_L1( xx,yy,n,plotopt )
%
% Performs a similar function to the fit_sinNtheta function, but minimises
% the error using the L1 norm. Solves by transformation to the linear
% programming problem
% This script finds a minimum length  fit to the data for a function of the
% form: y = A sin (w*x + phi)
% where A is the amplitude, w is the angular frequency, phi is offset
%
% INPUTS
%   xx = series of times/x-values/angles in degrees
%   yy = series of amplitudes/y-values
%   n  = number of full cycles in 360 degrees. 360/n is the wavelength
%   We = vector of weights for each of the (xx,yy) pairs = 1./sd^2
%   plotopt = whether to plot (1) or not (0)
% OUTPUTS
%   A = parameter 1 = amplitude of the sine wave
%   phi= parameter 2 = offset of the sine wave
%   rmserrorY = root mean square of the y errors
%   rmserrorX = attempted/rough root mean sqauare of x errors

% d = A*sin(w*x + phi)
% d = A*sin(phi)*cos(w*x) + A*cos(phi)*sin(w*x)
% d = m1*cos(w*t_i) + m2*sin(w*t_i)
% m1 = A*sin(phi)
% m2 = A*cos(phi)

% w=360/n;
if nargin<5, plotopt=0; end
if nargin<4, We=eye(length(xx)); end

% treat 4th argument as plotopt if it is just a 1 or zero
if nargin==4 && length(We)==1 && length(xx)>1
    plotopt=We; We=eye(length(xx));
elseif nargin>3
    We=diag(We);
end

M=2;
N=length(xx);
if length(yy)~=N, error('xx and yy must be the same length'); end
if size(xx,1)==1, xx=xx'; end
if size(yy,1)==1, yy=yy'; end

G=zeros(N,2);
for ix=1:N
G(ix,1) = cosd(n*xx(ix));
G(ix,2) = sind(n*xx(ix));
end

% L1 solution
% set up for linear programming problem
% min f*x subject to A x <= b,  Aeq x = beq

% variables
% m = mp - mpp
% x = [mp', mpp', alpha', x', xp']
% with mp, mpp length M and alpha, x, xp, length N
L = 2*M+3*N;
x = zeros(L,1);

% f is length L
% minimize sum aplha(i)/sd(i)
f = zeros(L,1);
f(2*M+1:2*M+N)=1.*sqrt(diag(We));

% equality constraints
Aeq = zeros(2*N,L);
beq = zeros(2*N,1);
% first equation G(mp-mpp)+x-alpha=d
Aeq(1:N,1:M)             =  G;
Aeq(1:N,M+1:2*M)         = -G;
Aeq(1:N,2*M+1:2*M+N)     = -eye(N,N);
Aeq(1:N,2*M+N+1:2*M+2*N) =  eye(N,N);
beq(1:N)                 =  yy;
% second equation G(mp-mpp)-xp+alpha=d
Aeq(N+1:2*N,1:M)               =  G;
Aeq(N+1:2*N,M+1:2*M)           = -G;
Aeq(N+1:2*N,2*M+1:2*M+N)       =  eye(N,N);
Aeq(N+1:2*N,2*M+2*N+1:2*M+3*N) = -eye(N,N);
beq(N+1:2*N)                   =  yy;

% inequality constraints A x <= b
% part 1: everything positive
A = zeros(L+2*M,L);
b = zeros(L+2*M,1);
A(1:L,:) = -eye(L,L);
b(1:L) = zeros(L,1);
% part 2; mp and mpp have an upper bound.  Note: 
% Without this constraint, a potential problem is that
% mp and mpp are individually unbounded, though their
% difference, m=mp-mpp, is not.  Thus there might be cases
% where the algroithm strays to very large mp and mpp.
A(L+1:L+2*M,:) = eye(2*M,L);
mls = (G'*G)\(G'*yy); %L2
mupperbound=10*max(abs(mls)); % might need to be adjusted for problem at hand
b(L+1:L+2*M) = mupperbound;

% solve linear programming problem
[x, fmin] = linprog(f,A,b,Aeq,beq);
fmin=-fmin;
mest = x(1:M) - x(M+1:2*M);
dpre = G*mest;

A = norm(mest);
phi = r2d(atan2(mest(1),mest(2)));
yyest=A*sind(n*xx + phi);
rmserrorY=norm((yy-yyest).*sqrt(diag(We)))/sqrt(N);

% L2 norm for comparison 
mest_L2 = (G'*We*G)\G'*We*yy; % Least squares comparison

A_L2 = norm(mest_L2);
phi_L2 = r2d(atan2(mest_L2(1),mest_L2(2)));
yyest_L2=A_L2*sind(n*xx + phi_L2);
rmserrorY_L2=norm((yy-yyest_L2).*sqrt(diag(We)))/sqrt(N);

% ugly stuff to get error in x
yy1=min(abs(yy),0.99999*A).*sign(yy);

xxest=zeros(N,2*n + 4);
xxest(:,1)=mod((asind(yy1/A)-phi),360)/n;
xxest(:,2)=mod(((180-asind(yy1/A))-phi),360)/n;
for ik=2:n
    xxest(:,2*ik - 1) = xxest(:,1) + (ik-1)*360/n;
    xxest(:,2*ik)     = xxest(:,2) + (ik-1)*360/n;
end
xxest(:,2*n + 1) = xxest(:,1) - 360/n;
xxest(:,2*n + 2) = xxest(:,2) - 360/n;
xxest(:,2*n + 3) = xxest(:,1) + 360/n + 360;
xxest(:,2*n + 4) = xxest(:,2) + 360/n + 360;

xxcomp=zeros(N,2*n + 4);
for ik=1:2*n + 4
    xxcomp(:,ik)=(xxest(:,ik)-xx).^2;
end
xxcomp=min(xxcomp,[],2);
rmserrorX=sqrt(mean(xxcomp));



% plot
if plotopt==1figure(99)
plot(xx,yy,'b.')
hold on
xxx=[0:1:360];
yyy=A*sind(n*xxx + phi);
yyy_L2=A_L2*sind(n*xxx + phi_L2);
plot(xxx,yyy,'r-')
plot(xxx,yyy_L2,'k--')
% plot(xxest,yy,'g.')
hold off
xlim([0 360])
title('blue=dobs red=L1 fit black=L2 fit')

end
