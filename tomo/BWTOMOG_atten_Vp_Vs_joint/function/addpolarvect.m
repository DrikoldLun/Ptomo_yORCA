function [ lat,lon,mag ] = addpolarvect( lat1,lon1,mag1,lat2,lon2,mag2 )
%[ lat,lon,mag ] = addpolarvect( lat1,lon1,mag1,lat2,lon2,mag2 )
% this function adds two vectors in polar coordinates
% input lat/lon must be in degrees
%
% Written by Zach Eilon, 2011

colat1=90-lat1;
colat2=90-lat2;

% Transform into cartesia n
x1=mag1*sind(colat1)*cosd(lon1);
y1=mag1*sind(colat1)*sind(lon1);
z1=mag1*cosd(colat1);
x2=mag2*sind(colat2)*cosd(lon2);
y2=mag2*sind(colat2)*sind(lon2);
z2=mag2*cosd(colat2);

x=x1+x2;
y=y1+y2;
z=z1+z2;

lon=r2d(atan2(y,x));
mag=sqrt(x^2 + y^2 + z^2);
colat=r2d(acos(z/mag));
% colat=r2d(atan2(sqrt(x^2 + y^2),z));
lat=90-colat;
end

