function [ A, phi, rmserrorY, rmserrorX ] = fit_sinNtheta( xx,yy,N,We,plotopt )
%function [ A, phi, rmserrorY, rmserrorX ] = fit_sinNtheta( xx,yy,N,We,plotopt )
%OR 
%function [ A, phi, rmserrorY, rmserrorX ] = fit_sinNtheta( xx,yy,N,plotopt )
%
% This script finds a least squares fit to the data for a function of the
% form: y = A sin (w*x + phi)
% where A is the amplitude, w is the angular frequency, phi is offset
%
% INPUTS
%   xx = series of times/x-values/angles in degrees
%   yy = series of amplitudes/y-values
%   N  = number of full cycles in 360 degrees. 360/N is the wavelength
%   We = vector of weights for each of the (xx,yy) pairs
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

% w=360/N;
if nargin<5, plotopt=0; end
if nargin<4, We=eye(length(xx)); end

% treat 4th argument as plotopt if it is just a 1 or zero
if nargin==4 && length(We)==1 && length(xx)>1
    plotopt=We; We=eye(length(xx));
elseif nargin>3
    We=diag(We);
end

nx=length(xx); ny=length(yy);
if nx~=ny, error('xx and yy must be the same length'); end

G=zeros(nx,2);

for ix=1:nx
G(ix,1) = cosd(N*xx(ix));
G(ix,2) = sind(N*xx(ix));
end
m = (G'*We*G)\G'*We*yy;
A = norm(m);
phi = r2d(atan2(m(1),m(2)));

yyest=A*sind(N*xx + phi);
rmserrorY=norm(yy-yyest)/sqrt(nx); % WAS

rmserrorY=norm((yy-yyest).*sqrt(diag(We)))/sqrt(nx);
rmserrorY=norm((yy-yyest).*sqrt(diag(We)/mean(diag(We))))/sqrt(nx);



% ugly stuff to get error in x
yy1=min(abs(yy),0.99999*A).*sign(yy);

xxest=zeros(nx,2*N + 4);
xxest(:,1)=mod((asind(yy1/A)-phi),360)/N;
xxest(:,2)=mod(((180-asind(yy1/A))-phi),360)/N;
for ik=2:N
    xxest(:,2*ik - 1) = xxest(:,1) + (ik-1)*360/N;
    xxest(:,2*ik)     = xxest(:,2) + (ik-1)*360/N;
end
xxest(:,2*N + 1) = xxest(:,1) - 360/N;
xxest(:,2*N + 2) = xxest(:,2) - 360/N;
xxest(:,2*N + 3) = xxest(:,1) + 360/N + 360;
xxest(:,2*N + 4) = xxest(:,2) + 360/N + 360;

xxcomp=zeros(nx,2*N + 4);
for ik=1:2*N + 4
    xxcomp(:,ik)=(xxest(:,ik)-xx).^2;
end
xxcomp=min(xxcomp,[],2);
rmserrorX=sqrt(mean(xxcomp));

% plot
if plotopt==1
figure(99)
plot(xx,yy,'b.')
hold on
xxx=[0:1:360];
yyy=A*sind(N*xxx + phi);
plot(xxx,yyy,'r-')
% plot(xxest,yy,'g.')
hold off
xlim([0 360])
end

