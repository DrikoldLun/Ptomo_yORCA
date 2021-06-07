function sol=pdeplate1(zdep0, tar0)
% function sol=pdeplate1(zdep, tar)
% solve the heat flow equation, 1D, for T=0 on top, T=1450 at Z=inf
% using the pdpde MATLAB solver.  
% Paramaters within are from GDH1, BUT allows for variable thermal
% conductivity - see pdex2pde to change
%
% INPUT:  zdep = depths to solve (km)
%         tar  = times to solve (Ma).   Both are arrays of solution
% OUTPUT:  
%      sol(I,J) gives T(tar(I),zdep(J))  in deg C
% GAA 2/13
ma2sec=3.156E+13;  % sec/Ma
if (nargin<2)
    zdep=0:1:150;
    tar=[0.1 0.2 0.4 1.0 2.0 4.0 8.0 16.0 32.0];
else
    zdep=zdep0;
    tar=tar0;
end
zdep=zdep.*1000;
tar=tar.*ma2sec;

sol=pdepe(0, @pdex2pde, @pdex2ic,@pdex2bc,zdep,tar);
sol = sol.*1421.5;    % correct to potential temperature for GDH1 (later add adiabat)
return

% --------------------------------------------------------------------------

function [c,f,s] = pdex2pde(x,t,u,DuDx)
% This is the main differential equation:  DuDx = DT/DZ  etc.
% C is rho*Cp
% F is K*DuDx    K=conductivity
% S = 0
%   GDH1 parameters
%  THIS VERSION allows K to increase by some factor in upper layer
cp=1.171E3;
kk=3.138;
rhom=3330;
c = cp*rhom;
zmin=5000;
kfac = 5.0;    % if depth<zmin then K= kk*kfac

f=kk*DuDx;

% Jack up conductivity in upper layer  (zmin meters) by factor kfac
% Next could jack up conductivity for shallow
if x <= zmin, f = f*kfac;  end
s=0;
return

% --------------------------------------------------------------------------

function u0 = pdex2ic(x)
% Set the Initial Condition  - to 1450C = 1
if x==0
  u0 = 0;
else
  u0 = 1;
end
return

% --------------------------------------------------------------------------

function [pl,ql,pr,qr] = pdex2bc(xl,ul,xr,ur,t)
% Set the boundary condition  on left, right side  (FIX THEM)
pl = ul;  % pdepe will use pl=ql=0 because this 
ql = 0;  % problem has a singularity: xl=0 and m>0
pr = ur-1;
qr = 0;    % hopefully this sets zero flux at infinity?
return