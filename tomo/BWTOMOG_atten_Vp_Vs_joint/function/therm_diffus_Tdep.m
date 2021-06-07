function [kappa,k,Cp,rho]  = therm_diffus_Tdep( T,rho )
%k  = therm_conduct_Tdep( T )
%   
% Function to calculate the temperature-dependent thermal diffusivity,
% accounting for thermal dependence of conductivity (by both radiative and phonon heat transfer), heat capacity, and density. All equations are from
% McKenzie et al. 2005 EPSL
% 
% INPUTS:
%   T       = temperature, in deg. centigrade
%   rho     = density, in kg/m^3 [optional, default is to calc. it from temp and rho0=3300]
% 
% OUTPUTS:
%   kappa   = diffusivity, m^2 / s
%   k       = conductivity,  W / m / K
%   Cp      = heat capacity, J / kg / K
%   rho     = density, kg / m^3
% 
%  Recall:  kappa = k/(rho*Cp)
% 
%  ZCE 07/2016

Tk = T+273;

%% conductivity
% for Olivine, can use expression of Xu et al (McKenzie equation 6)
k = 4.08.*(298/Tk).^0.406;

%% heat capacity
% (McKenzie equation 10)
Cp = 233.18 - 1801.6.*Tk.^(-0.5) - 26.794e7.*Tk.^-3;

%% density
if nargin<2 || isempty(rho)
    rho = 3330 * exp( - (2.832e-5.*T + 0.5*3.79e-8*(Tk.*Tk - 273*273) ) );
end

kappa = k./(rho.*Cp);

end

