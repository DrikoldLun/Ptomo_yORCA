function [foraz,backaz,gcdist,gcarc] = ellipsep (lat1,long1,lat2,long2)
% [foraz,backaz,gcdist,gcarc] = ellipsep (lat1,long1,lat2,long2)
% This function explicitly  computes the geographical distance and forward
% and back azimuths between point 1 and point 2 on an ellipsoidal Earth,
% using Vincenty's formulae and iterating. In theory, it is accurate to
% within 0.5mm. see http://en.wikipedia.org/wiki/Vincenty%27s_formulae
% 
% Useful for calculating raypaths from events (point 1) to stations (point
% 2).
% 
% Inputs: coordinates (in degrees) of point 1 (lat1,long1) and point 2
% (lat2,long2) 
% Outputs: forward azimuth (foraz) from 1==>2, backazimuth (backaz) from
% 2==>1, great circle distance (gcdist) - distance over the surface from 1
% to 2, great circle arc (gcarc) - angle bisected at the centre of the
% earth, angular separation. 
% 
% 
% a     length of major axis of the ellipsoid (radius at equator);	(6378137.0 metres in WGS-84)
% f     flattening of the ellipsoid;	(1/298.257223563 in WGS-84)
% b = (1 - f)*a     length of minor axis of the ellipsoid (radius at the poles);
% ?1, ?2 = lat1,lat2	latitude of the points;
% U1 = arctan[(1 ? ?) tan ?1],
% U2 = arctan[(1 ? ?) tan ?2]	reduced latitude (latitude on the auxiliary sphere)
% L = L2 - L1	difference in longitude of two points;
% ?1, ?2 = long1,long2	longitude of the points on the auxiliary sphere;
% ?1, ?2 = az1,az2      forward azimuths at the points;
% ?	= alpha     azimuth at the equator;
% s	= gcdist    ellipsoidal distance between the two points;
% ?	= gcarc     arc length between points on the auxiliary sphere;
% 
% Written by Zach Eilon, 2012

if lat1>lat2
    sign=-1;
else
    sign=1;
end

d2r=pi/180;
lat1=d2r*lat1; %degrees2radians
long1=d2r*long1; %degrees2radians
% long1=mod(long1,2*pi); % convert to 0-long1-2pi
lat2=d2r*lat2; %degrees2radians
long2=d2r*long2; %degrees2radians
% long2=mod(long2,2*pi); % convert to 0-long2-2pi

a=6378137.0;
f=(1/298.257223563);
b=(1-f)*a;
U1=atan((1-f)*tan(lat1)); % Reduced latitude
U2=atan((1-f)*tan(lat2));

L=abs(long1-long2); % L is absolute distance, as the longs have been converted to 0?long?2pi
% if L > pi 
%     L = 2*pi - L; % so 0?L?pi
% end

% L=long1-long2
lambda=L;
Dlambda=1;
while Dlambda > 1e-12
sin_s = sqrt((cos(U2)*sin(lambda))^2 + (cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(lambda))^2);
cos_s = sin(U1)*sin(U2) + cos(U1)*cos(U2)*cos(lambda);
s = atan2(sin_s,cos_s); % N.B. here 's' is sigma=gcarc - don't confuse with ellipsoidal distance=gcdist, above. 
sin_a = cos(U1)*cos(U2)*sin(lambda)/sin_s;
cos_a = sqrt(1 - sin_a^2);
cos2_a= cos_a^2;
cos_2sm = cos_s - 2*sin(U1)*sin(U2)/cos2_a;
C=(f/16)*cos2_a*(4 + f*(4 - 3*cos2_a));
lambdanew = L + (1 - C)*f*sin_a*(s + C*sin_s*(cos_2sm + C*cos_s*(-1 + 2*cos_2sm^2)));
Dlambda=abs(lambda-lambdanew);
lambda=lambdanew;
end
u2 = cos2_a*(a^2 - b^2)/b^2;
A = 1 + (u2/16384)*(4096 + u2*(-768 + u2*(320 - 175*u2)));
B = (u2/1024)*(256 + u2*(-128 + u2*(74 - 47*u2)));
Ds = B*sin_s*(cos_2sm + (B/4)*(cos_s*(-1 + 2*cos_2sm^2) - (B/6)*cos_2sm*(-3 + 4*sin_s^2)*(-3 + 4*cos_2sm^2)));
gcarc = s/d2r; % in degrees
gcdist = b*A*(s - Ds); % in metres
foraz = (1/d2r)*atan2(cos(U2)*sin(lambda),(cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(lambda))); %degrees from N
backaz = mod(180+(1/d2r)*atan2(cos(U1)*sin(lambda),(-sin(U1)*cos(U2) + cos(U1)*sin(U2)*cos(lambda))),360); %degrees from N
if backaz > 180; backaz=backaz-360;
end