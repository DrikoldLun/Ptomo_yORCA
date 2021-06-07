function [ x, y ] = lonlat2xy( olon, olat, lon, lat )
% [ x, y ] = lonlat2xy( olon, olat, lon, lat )
% 
% Convert lon/lat (deg) to x/y (m) using the reference ellipsoid

global ellipsoidGRS80
if isempty(ellipsoidGRS80)
    ellipsoidGRS80 = referenceEllipsoid('wgs84');
end

[x, y, ~] = geodetic2enu(lat, lon, zeros(size(lon)), olat, olon, 0, ellipsoidGRS80);

% If no mapping toolbox, use this version
% x = -(olon - lon)*cos(olat*pi/180)*111.1 * 1000; % m
% y = -(olat - lat)*111.1 * 1000; % ms



end

