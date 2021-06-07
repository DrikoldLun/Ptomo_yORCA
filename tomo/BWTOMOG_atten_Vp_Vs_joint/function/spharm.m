function [Ylm] = spharm(L,M,THETA,PHI,cs,plotopt)
% [Ylm] = spharm(L,M,cs,THETA,PHI,plotopt)
% function to calculate the value of a spherical harmonic at a given point
% on the Earth
%
% INPUTS
% L - the angular degree
% M - the azimuthal order
% THETA - a vector of LATitudes, in degrees from 90 (N) to -90 (S)
% PHI - a vector of LONgitudes, in degrees from 180 (E) to -180 (W)
% cs - the odd/even mode - 'c' is cos, real; 's' is sin, imag.
% plotopt - option to plot on the spherical harmonic patches. 
%
% OUTPUTS
% Ylm - a vector of amplitude of the function at the points described by
% theta and phi. By default this is complex.

if M>L
    error('The degree (l) must be greater than or equal to the order (m)')
end

if nargin > 4 && isnumeric(cs)==1 % skipped c/s option
	plotopt = cs;
end
THETA = d2r(ones(size(THETA))*90 - THETA);
PHI = d2r(mod(PHI,360));

P = legendre(L,cos(THETA));
a1 = (2*L + 1)/(4*pi);
a2 = factorial(L-M)/factorial(L+M);
a3 = abs((-1)^M)*sqrt(a1*a2);
E = exp(1i*M*PHI);
if nargin > 4 && isnumeric(cs)~=1
if strcmp(cs,'c')
    E = real(E);
elseif strcmp(cs,'s')
    E = imag(E);
end
end
Ylm = a3*P(M+1)*E;
