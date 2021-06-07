function [ A1, phi1, A2, phi2, rmserrorY, rmserrorX ] = fit_sinNtheta( xx,yy,N,We,plotopt )
%function [ A, phi, rmserrorY, rmserrorX ] = fit_sinNtheta( xx,yy,N,We,plotopt )
%OR 
%function [ A, phi, rmserrorY, rmserrorX ] = fit_sinNtheta( xx,yy,N,plotopt )
%
% This script uses the least squares method (like in Chevrot 2000) to find
% a fit for a function of the form:
% y = A1 sin (x + phi1) + A2 sin (2x + phi2) 
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

% d = A1*sin(x + phi1) + A2*sin(2x + phi2)
% d = A1*sin(phi1)*cos(x) + 1*cos(phi1)*sin(x) + A2*sin(phi2)*cos(2x) + A2*cos(phi2)*sin(2x)
% d = m1*cos(x) + m2*sin(x) + m3*cos(2x) + m4*sin(2x)
% m1 = A1*sin(phi1)     A1 = ?(m1^2 + m2^2)
% m2 = A1*cos(phi1)     phi1 = atan(m1/m2)
% m3 = A2*sin(phi2)     A2 = ?(m3^2 + m4^2)
% m4 = A2*cos(phi2)     phi2 = atan(m3/m4)

% construct equations   Gm = d  where
% G(i,:) = [ cos(x_i) sin(x_i) cos(2*x_i) sin(2*x_i) ]
% m = [ m1; m2; m3; m4 ]
% d = dobs (column vector)
% and solve using least squares

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

G=zeros(nx,4);

for ix=1:nx
G(ix,1) = cosd(xx(ix));
G(ix,2) = sind(xx(ix));
G(ix,3) = cosd(2*xx(ix));
G(ix,4) = sind(2*xx(ix));
end
mest = (G'*We*G)\G'*We*yy;
A1 = norm(mest(1:2));
A2 = norm(mest(3:4));
phi1 = r2d(atan2(mest(1),mest(2)));
phi2 = r2d(atan2(mest(3),mest(4)));

yyest = A1*sind(xx + phi1) + A2*sind(2*xx + phi2);
rmserrorY=norm(yy-yyest)/sqrt(nx); % WAS
rmserrorY=norm((yy-yyest).*sqrt(diag(We)))/sqrt(nx);
rmserrorY=norm((yy-yyest).*sqrt(diag(We)/mean(diag(We))))/sqrt(nx);


% % ugly stuff to get error in x
% yy1=min(abs(yy),0.99999*A).*sign(yy);
% 
% xxest=zeros(nx,2*N + 4);
% xxest(:,1)=mod((asind(yy1/A)-phi),360)/N;
% xxest(:,2)=mod(((180-asind(yy1/A))-phi),360)/N;
% for ik=2:N
%     xxest(:,2*ik - 1) = xxest(:,1) + (ik-1)*360/N;
%     xxest(:,2*ik)     = xxest(:,2) + (ik-1)*360/N;
% end
% xxest(:,2*N + 1) = xxest(:,1) - 360/N;
% xxest(:,2*N + 2) = xxest(:,2) - 360/N;
% xxest(:,2*N + 3) = xxest(:,1) + 360/N + 360;
% xxest(:,2*N + 4) = xxest(:,2) + 360/N + 360;
% 
% xxcomp=zeros(nx,2*N + 4);
% for ik=1:2*N + 4
%     xxcomp(:,ik)=(xxest(:,ik)-xx).^2;
% end
% xxcomp=min(xxcomp,[],2);
% rmserrorX=sqrt(mean(xxcomp));

% plot
if plotopt==1
figure(99)
plot(xx,yy,'b.')
hold on
xxx=[0:1:360];
yyy=A1*sind(xxx + phi1) + A2*sind(2*xxx + phi2);
plot(xxx,yyy,'r-')
% plot(xxest,yy,'g.')
hold off
xlim([0 360])
end

