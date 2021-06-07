function [ P ] = P_from_Z_rho( Z,rho )
% [ P ] = P_from_Z_rho( Z,rho )
%  Calculate pressure from depth and density vectors (or a density scalar);
% Z    in km
% rho  in kg/m^3
% P    in Pa 
% N.B. the integration requires the spacing of Z to be small with respect
% to the radius of the Earth such that g does not change much over the
% pressure increments.


Z = Z(:); rho = rho(:); % make columns
if isscalar(rho), rho = rho*ones(size(Z)); end

Me = 5.972e24;
Re = 6371e3;
G = 6.67e-11;

if Z(1)~=0
    Z = [0;Z];
    rho = [3300;rho];
end
Z = Z*1000; % in m
R = Re-Z;
P = zeros(size(Z)); % pressure of each point
M = zeros(size(Z)); % mass inside that radius

% start at outside, R=Re, M=Me
M(1) = Me;
for ip = 2:length(Z)
    dz = Z(ip)-Z(ip-1);
    avrho = mean(rho([ip-1:ip]));
    avR = mean(R([ip-1:ip]));
    dM = 4*pi*avR^2*dz*avrho;
    M(ip) = M(ip-1) - dM;
    avg = G * mean(M([ip-1:ip])) * avR.^-2; % acc. of gravity
    dP = avrho * dz * avg;
    P(ip) = P(ip-1) + dP;
    
end


