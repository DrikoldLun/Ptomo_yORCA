function [ datZRT ] = psvsh2zrt( datPSVSH, a, B, p, zsgn )
% [ datZRT ] = psvsh2zrt( datPSVSH, a, B, phi, zsgn )
%   Function to transform P-SV-SH into Z-R-T  seismogram components using
%   free-surface transfer matrix (Kennett, 1991) from Abt et al 2010
% INPUTS:
% datPSVSH  - a Nx3 matrix of particle motions in the order P, SV, SH (SV +ve outwards)
% a         - near surface compressional wavespeed (km/s)
% B         - near surface shear wavespeed (km/s)
% p         - ray parameter = sin(inc)/v where v is a or B   (in s/km)
% zsgn      - Z comp sign convention (+1 if Z +ve up, -1 if Z +ve down)

if nargin<5 || isempty(zsgn)
    zsgn = 1;
end

if a*p > 1, error('P-conversion is post-critical!'); end

xa = sqrt(1 - (a*p)^2);
xB = sqrt(1 - (B*p)^2);

M = [-zsgn*(B^2 * p^2 - 0.5)/xa      ((p*B^2)/a)      0
            -zsgn*p*B          (0.5 - B^2 * p^2)/xB   0
             0                      0          0.5];

invM = inv(M);

datZRT = datPSVSH*invM';

end

