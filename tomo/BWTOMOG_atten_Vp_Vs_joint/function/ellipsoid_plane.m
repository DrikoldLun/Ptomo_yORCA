% function [ v1,v2,faz_out ] = ellipsoid_plane( faz,ap,bp,va,vb,vc,Z )
%[ v1,v2,faz_out ] = ellipsoid_plane( faz,ap,bp,va,vb,vc,Z )
% NB STILL PROBLEMS WITH THIS FUNCTION!!!
%   function to calculate the minor and major axes of the ellipse produced
%   from the intersection of a plane (with normal Z) and an ellipsoid with
%   arbitrary orientation.
% 
%   Ellipsoid shape is given by a, b, c, the orthogonal lengths of the
%   maximum and minimum axes, in clockwise orientation. Ellipsoid
%   orientation is given by a set of Euler angles, which rotate it with
%   respect to the coordinate system E,N,Z (Z is up, so clockwise system)
% 
%  Following all conventions from Gendzwill and Stauffer, 1981.
% 
% INPUTS:
%   faz     azimuth of a axis
%   ap      plunge of a axis
%   bp      plunge of b axis (not true plunge from vertical, rotation about a axis
%   va      major axis of ellipsoid
%   vb      intermediate axis of ellipsoid
%   vc      minor axis of ellipsoid
%   Z       [e,n,z] vector of normal to a plane in the e,n,z coord system
% 
%  N.B. require va >= vb >= vc > 0
% 
% OUTPUTS:
%   v1      major axis of plane section ellipse 
%   v2      minor axis of plane section ellipse
%   faz_out rake of major axis from Z foraz (clk-wise from foraz to fast)

%% FOR TESTING
faz = 0;
ap = 1;
bp = 0;
va = 2;
vb = 1.5;
vc = 1;
Z = [0,tand(1),1];
%% for testing end

Z = Z/norm(Z);



% have to translate faz,ap,bp to Euler angles (which are rotations from
% different reference axes... (i.e. E not N)
phi = 90-faz;% 1st Euler angle
the = ap;
psi = bp;

J = [ cosd(the)*cosd(phi)*cosd(psi) - sind(phi)*sind(psi)
      cosd(the)*sind(phi)*cosd(psi) + cosd(phi)*sind(psi)
     -sind(the)*cosd(psi) ]
K = [-cosd(the)*cosd(phi)*sind(psi) - sind(phi)*cosd(psi)
     -cosd(the)*sind(phi)*sind(psi) + cosd(phi)*cosd(psi)
      sind(the)*sind(psi) ]
L = [ sind(the)*cosd(phi)
      sind(the)*sind(phi)
      cosd(the) ]
  
aa = acosd(  dot(L,Z)/(norm(L)*norm(Z))  ); % alpha in paper

if dot(cross(L,Z),K)==0 % gotta account for atan2(0,0)=0 (we want 90)
    bb = 90;
else
    bb  = r2d(atan2(dot(cross(cross(L,Z),K),Z),dot(cross(L,Z),K))); % beta in paper
end

ca = cosd(aa); cca = cosd(aa)^2;
sa = sind(aa); ssa = sind(aa)^2;
cb = cosd(bb); ccb = cosd(bb)^2;
sb = sind(bb); ssb = sind(bb)^2;

F1 = (cca*ccb)/(va^2) + (cca*ssb)/(vb^2) + ssa/(vc^2);
F2 = 2*ca*cb*sb*(vb^-2 - va^-2);
F3 = ssb/(va^2) + ccb/(vb^2);

        % this is given as plus-minus in the paper...
if F2==0 %% should this have a condition for F3-F1 too?
    tang0 = [Inf -Inf]';
else
    tang0 = [(F3-F1)/F2 + sqrt(((F3-F1)/F2)^2 + 1)
             (F3-F1)/F2 - sqrt(((F3-F1)/F2)^2 + 1)];
end
g0 = atand(tang0); % gamma in paper
% g0 = g0(1); % choose plus...?

P = ( F1.*cosd(g0).^2 + F2.*cosd(g0).*sind(g0) + F3.*sind(g0).^2 ).^-0.5;
S = ( F3.*cosd(g0).^2 - F2.*cosd(g0).*sind(g0) + F1.*sind(g0).^2 ).^-0.5;

% NOTE: These appear to be in the wrong order in the paper, i.e. here S>P
v1 = S;
v2 = P;

% H tied only to azimuth of Z - ref we neeed for faz_out
% this is the horizontal vector in the plane of the section, 90 deg
% anticlockwise from foraz
H = cross([0,0,1],Z) 

% Y,X tied to orie of indicatrix, and azimuth of Z:
Y = cross(L,Z)
X = cross(Y,Z) % slopes 'down' the section plane

ry = r2d(atan2(dot(cross(Y,H),Z),dot(Y,H)));

rp = 90 - g0 + ry % Clockwise angle within plane from H to v1 dir.
% we want the angle from foraz to v1 dir, = rp-90
faz_out = rp-90;

faz_out = mod(faz_out+90,180)-90; %transform to -90 to +90

% Pick correct orie - the one for which v1>=v2
[~,ab] = max(v1./v2);
v1 = v1(ab)
v2 = v2(ab)
faz_out = faz_out(ab)





% end

