function [lat2, lon2] = unreckon(lat1, lon1, gcarc, backaz)
% [latout, lonout] = unreckon(lat, lon, gcarc, backaz)
% For a given location, this function calculates the position on the globe
% that will have the desired backazimuth to the original location
% all angles and arc distances in degrees
acc=0.0001;
theta1=90-lat1;

dlon=asind(sind(gcarc)*sind(backaz)/sind(theta1));
lon2=lon1-dlon;


cG  = cosd(gcarc);
sGcB= sind(gcarc)*cosd(backaz);
cT  = cosd(theta1);
sTcP= sind(theta1)*cosd(dlon);

A = cG^2 - cT^2;
B = cG*sTcP - cT*sGcB;
% sind(theta2)=A/B
C = cG*sGcB - cT*sTcP;
D = cT*sGcB - cG*sTcP;
% cosd(theta2)=C/D
theta2=atan2(A*D,B*C);
theta2=r2d(theta2);
theta2=mod(theta2,180);
lat2=90-theta2;

%check
[rng,az]=distance(lat2,lon2,lat1,lon1);
if (rng-gcarc)>acc || (az-backaz)>acc
    error('PROBLEM: no solution for these parameters')
end

end

