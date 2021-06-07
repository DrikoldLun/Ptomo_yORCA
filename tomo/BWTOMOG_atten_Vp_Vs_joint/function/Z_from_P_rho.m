function [ Z ] = Z_from_P_rho( P,rho )
% [ Z ] = Z_from_P_rho( P,rho )
%  Calculate depth from pressure and density vectors (or a density scalar);
% P    in Pa
% rho  in kg/m^3
% Z    in km
% N.B. the integration requires the spacing of P to be small with respect
% to the pressure at the centre of the Earth such that g does not change
% much over the pressure increments.


P = P(:); rho = rho(:); % make columns
if isscalar(rho), rho = rho*ones(size(P)); end

Me = 5.972e24;
Re = 6371e3;
G = 6.67e-11;

if P(1)~=0
    P = [0;P];
    rho = [3300;rho];
end

Z = zeros(size(P)); % depth of each point
R = zeros(size(P)); % radius of each point
M = zeros(size(P)); % mass inside that radius

% start at outside, R=Re, M=Me
R(1) = Re;
M(1) = Me;
for ip = 2:length(P)
    dP = P(ip)-P(ip-1);
    avrho = mean(rho([ip-1:ip]));
    g = G*M(ip-1)*R(ip-1).^-2; % acc. of gravity
    
    dg = 1e5;
    while abs(dg) >0.001; % iterate on gravity until it's representative of the middle of the layer
    dz = dP./(avrho*g);
    
    Z(ip) = Z(ip-1)+dz;
    R(ip) = R(ip-1)-dz;
    
    avR = R(ip) - dz/2;
    dM = 4*pi*avR^2*dz*avrho; % mass of the shell
    
    M(ip) = M(ip-1)-dM;
    g_ = G*mean(M([ip-1:ip]))*avR.^-2;
    dg = g-g_;
    g = g_;
    
    end
    
    
end

Z = Z/1000; % change from m to km

