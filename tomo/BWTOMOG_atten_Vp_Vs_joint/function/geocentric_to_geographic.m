function [olat] = geocentric_to_geographic(ilat,direction)
% [olat] = geocentric_to_geographic(ilat,direction)
% 
% Function to transform back and forth between geocentric latitude (angle
% subtended by line from point to centre of Earth) and the geographic
% latitude (angle made with the equatorial plane from a line tangent to the
% surface at a point - this is the one we conventionally plot).
% 
%  INPUTS:
%    ilat      - input latitude (in degrees)
%    direction - 'forward'=geoc2geog(default) or 'backward'=geog2geoc
% 
%  OUTPUT:
%    olat      - output latitude
% 
% If phic is the geocentric latitude and phig is the geographic latitude. 
%   tan(phig) * (1-f)^2 = tan(phic)     , where f is the ellipticity
% 
% Example: geocentric_to_geographic(40,'forward') transforms geocentric
% latitude of 40 to a geographic latitude of 40.1896.

if nargin < 2 || isempty(direction)
    direction = 'forward'; 
end


f = 1./298.25722; % Ellipticity of Earth, World Geodetic Survey (WGS-84)

fac = (1-f).^2;

switch direction
    case 'forward'
        olat = atand( tand(ilat)./fac );
    case 'backward'
        olat = atand( tand(ilat).*fac );
end