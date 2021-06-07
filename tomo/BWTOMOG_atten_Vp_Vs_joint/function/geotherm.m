function [ TempZ ] = geotherm( age,coolopt,Z,Tpot,kappa,dTdz_ad,Zplate  )
% [ TempZ ] = geotherm( age,coolopt,Z,Tpot,kappa,dTdz_ad,Zplate  )
%   function to calculate an oceanic geotherm based on a given cooling
%   model, based on Turcotte and Schubert, equations 4-125 and 4-130
% 
% INPUTS
%    age     = age of plate (Ma) = time of cooling = distance/spread_rate
% 
%    coolopt = cooling model option, 'plate' or 'halfspace' 
%               [optional, default is 'halfspace']
%    Z       = vector of depths (km) 
%               [optional, default is 0:5:150]
%    Tpot    = mantle potential temperature (C) 
%               [optional, default is 1300C]
%    kappa   = thermal diffusivity (m^2/s)      
%               [optional, default is 1e-6 m^2/s]
%               NB optionally include thermal dependence of diffusivity
%    dTdz_ad = adiabatic temperature gradient (degrees C/km)
%               [optional, default is 0.3 ?C/km]
%    Zplate  = thickness of plate if using plate model (km)
%               [optional, default is 125km]
% 
% OUTPUTS
%    TempZ   = vector of temperatures (in C) at each depth
% 
% Z. Eilon Jan 2015


if nargin < 2
    coolopt = 'halfspace';
end
if nargin < 3
    Z = [0:5:150]';
end
if nargin < 4
    Tpot = 1300;
end
if nargin < 5
    kappa = 1e-6;
end
if nargin < 6 
    dTdz_ad = 0.3;
end
if nargin < 7
    Zplate = 125;
end

% T0 = 273;
% T1 = Tpot + T0;

age = 365*24*3600*age*1e6; % time in s
Z = Z*1000; % depth in m
Zplate = Zplate*1000;

%% Adiabat
Tad = Tpot + Z*dTdz_ad/1000;

%% Cooling
if strcmp('halfspace',coolopt)

% TempZ = T0 + (T1-T0)*erf(Z/(2*sqrt(kappa*t))); % T&C eqn. 4-125
% shakes out as
TempZ = Tad.*erf(Z./(2*sqrt(kappa.*age)));

elseif strcmp('plate',coolopt)

%  T = T0 +  (T1-T0)*[(Z/Zplate) + (2/pi)SUM[(1/n)*exp()*sin()]
% shakes out as
T = (Z/Zplate).*ones(size(Z));
niter = 1e5;
for ii = 1:niter
    T = T + (2/pi)*(1/ii)*exp((-kappa*ii^2*pi^2*age)/(Zplate^2))*sin(ii*pi*Z/Zplate);
end
T(Z>=Zplate) = 1;
TempZ = Tad.*T;    
    
end

end

