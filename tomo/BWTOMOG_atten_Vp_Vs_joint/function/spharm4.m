function [Ymn,PHI,THETA,Xm,Ym,Zm]=spharm4(L,M,RES,PLOT_FLAG)
%This function generates the Spherical Harmonics basis functions of degree
% L and order M.
%
% SYNTAX: [Ymn,PHI,THETA,X,Y,Z]=spharm4(L,M,RES,PLOT_FLAG);
%
% INPUTS:
%
% L         - Spherical harmonic degree, [1x1]
% M         - Spherical harmonic order,  [1x1]
% RES       - Vector of # of points to use [#Phi x #Theta points],[1x2] or [2x1] 
% PLOT_FLAG - Binary flag to generates a figure of the spherical harmonic surfaces (DEFAULT=1)
%
%
% OUTPUTS:
% 
% Ymn   - Spherical harmonics coordinates, [RES(1) x RES(2)]
% PHI - Circumferential coordinates, [RES(1) x RES(2)]
% THETA   - Latitudinal coordinates, [RES(1) x RES(2)]
% X,Y,Z - Cartesian coordinates of magnitude, squared, spherical harmonic surface points, [RES(1) x RES(2)]
%
%
% NOTE: It is very important to keep the various usages of PHI and THETA
% straight.  For this function PHI is the Azimuthal/Longitude/Circumferential 
% coordinate and is defined on the interval [0,2*pi], whereas THETA is the 
% Altitude/Latitude/Elevation and is defined on the interval [0,pi].  Also note that 
% the conversion to cartesian coordinates requires that THETA be offset by pi/2 so 
% that the conversion is on the interval [-pi/2,pi/2].
%
% DBE 2005/09/30
% ZJE 2013/01/30 edited

% Define constants (REQUIRED THAT L(DEGREE)>=M(ORDER))
if nargin==0
  L=3;   % DEGREE
  M=2;   % ORDER
  RES=[55 55];
end

if nargin<3
  RES=[25 25];  
  PLOT_FLAG=1;
end

if nargin<4
  PLOT_FLAG=1;
end

if L<M, error('The ORDER (M) must be less than or eqaul to the DEGREE(L).'); end

PHI=linspace(0,2*pi,RES(1));  % Azimuthal/Longitude/Circumferential
THETA  =linspace(0,  pi,RES(2));  % Altitude /Latitude /Elevation
[PHI,THETA]=meshgrid(PHI,THETA);

Lmn=legendre(L,cos(THETA));

if L~=0
  Lmn=squeeze(Lmn(M+1,:,:));
end

a1=((2*L+1)/(4*pi));
a2=factorial(L-M)/factorial(L+M);
C=sqrt(a1*a2);

Ymn=(-1^M)*C*Lmn.*exp(1i*M*PHI);

[Xm,Ym,Zm]=sph2cart(PHI,THETA-pi/2,abs(Ymn).^2);
[Xr,Yr,Zr]=sph2cart(PHI,THETA-pi/2,real(Ymn).^2);
[Xi,Yi,Zi]=sph2cart(PHI,THETA-pi/2,imag(Ymn).^2);
% [Xp,Yp,Zp]=sph2cart(PHI,THETA-pi/2,angle(Ymn).^2);

if PLOT_FLAG
f=figure; axis off; hold on;
  axes('position',[0.0500 0 0.2666 1]); 
    surf(Xm,Ym,Zm); 
    axis equal off; %rot3d;
    light; lighting phong; camzoom(1.3);
  axes('position',[0.3666 0 0.2666 1]); 
    surf(Xr,Yr,Zr); 
    axis equal off; %rot3d;
    light; lighting phong; camzoom(1.3);
  axes('position',[0.6833 0 0.2666 1]); 
    surf(Xi,Yi,Zi); 
    axis equal off; %rot3d;
    light; lighting phong; camzoom(1.3);
  axes('position',[0 0.9 1 0.1]); axis off;
    t(1)=text(0.50,0.25,'Spherical Harmonics','HorizontalAlignment','Center');
  axes('position',[0 0 1 0.1]); axis off;
    t(2)=text(0.20,0.25,['|Y^',num2str(M),'_',num2str(L),'|^2'],'HorizontalAlignment','Center');
    t(3)=text(0.50,0.25,['Real(Y^',num2str(M),'_',num2str(L),')^2'],'HorizontalAlignment','Center');
    t(4)=text(0.80,0.25,['Imag(Y^',num2str(M),'_',num2str(L),')^2'],'HorizontalAlignment','Center');
%   subfig(gcf,10,5,12);
end

end