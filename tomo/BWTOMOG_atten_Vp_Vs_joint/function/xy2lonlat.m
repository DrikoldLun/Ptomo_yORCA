function [ lon, lat ] = xy2lonlat( olon, olat, x, y,az )
%[ lon, lat ] = xy2lonlat( olon, olat, x, y,az )
%
% Convert x/y (m) to lat/lon (degrees) using the reference ellipsoid
% option to rotate by azimuth "az" (defined as clockwise from N)

global ellipsoidGRS80

if nargin < 5 
    az = 0;
end

if az~=0
    x0 = x;
    y0 = y;
    x = x0.*cosd(az) - y0.*sind(az);
    y = x0.*sind(az) + y0.*cosd(az);
end

if isempty(ellipsoidGRS80)
    ellipsoidGRS80 = referenceEllipsoid('wgs84');
end

[lat, lon, ~] = enu2geodetic(x, y, zeros(size(x)), olat, olon, 0, ellipsoidGRS80);

end

