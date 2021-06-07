function [fastaz, DinvV, Vav, anis] = aniso3D(inth,raz,faz,ap,bp,anistype,aniscustom)
% [fastaz, DinvV, Vav, anis] = aniso3D(inth,raz,faz,ap,bp,anistype)
%
% function to output fast directions and velocity differences given a
% propagation vector and fabric orientation and type:
% 
% INPUTS
% Defining ray vector
%   inth    =  incidence angle 
%   raz     =  propagation azimuth 
% Defining fabric
%   faz     =  azimuth of fast axis
%   ap      =  plunge of the a-axis
%   bp      =  plunge of the b-axis 
%   anistype=  fabric type - see "anisotens.m" for details... 12 is ortho
%  TO DO CUSTOM ANISOTROPY: USE ANISTYPE==99 AND INCLUDE THE BELOW
%   aniscustom = [SVav,SVanis]: nominal Vs_average and Vs_anisotropy (%) in HEXAGOAL
% Note all angles are in degrees
% 
% OUTPUTS
%   fastaz  = azimuth (cw from N) of surface projection of fast direction
%   DinvV   = inverse shear wave velocity diff (s/km) (1/V_slow - 1/V_fast) 
%   Vav     = average shear wave velocity 0.5*(V_slow + V_fast)
%   anis    = anisotropy strength, as a percentage (V_fast-V_slow)/(V_fast+V_slow)
%  
% Define some directions
% (1) is global vertical    or      fabric a-axis
% (2) is global east        or      fabric b-axis
% (3) is global north       or      fabric c-axis
%
% Written by Zach Eilon, 2013


%% Defining ray vector
% % NOT clockwise difference between fast and propagation angle - i.e. if fast is zero and foraz =10, az=10.
rvec = zeros(3,1);
rvec(1) = cosd(inth);           % Vert (up is +ive)
rvec(2) = sind(inth)*sind(raz);  % East
rvec(3) = sind(inth)*cosd(raz);  % North
% rvec is propagation vector in coords (1)=fast (2)=90?cw-horiz from fast (3)=90?cw-vert from fast 
%% Defining fabric direction
rot0 = [0 0 1; 0 1 0; 1 0 0]; % to move (a) axis to (3) direction
rot1 = [1  0  0 ; 0 cosd(faz) -sind(faz) ; 0 sind(faz) cosd(faz)];
rot2 = [cosd(ap) 0 sind(ap); 0  1  0; -sind(ap) 0 cosd(ap)];
rot3 = [cosd(-bp)  -sind(-bp)  0 ; sind(-bp) cosd(-bp) 0 ; 0  0  1]; % "-" because anti-clockwise to plunge b
rotall = rot0*rot1*rot2*rot3;
% rotall is rotation from global (1)(2)(3) to local (a)(b)(c)

if anistype==99
    SVav = aniscustom(1); SVanis=aniscustom(2);
    % based on from Anderson 1989 "Transverse isotropic" (2)  11.6%S
    A=205.2; C=314.1; F=69.6; L=72.3; N=62.3; 
    den=3.3;
    N = SVav^2 * den * (1 - 0.01*SVanis)^2;
    L = SVav^2 * den * (1 + 0.01*SVanis)^2;
    [ cc,~ ] = Etensor_hex( A,C,F,L,N,1 );
else
    [cc, CC6, den] = anisotens(anistype);
end

n=rotall*rvec; % propagation vector in fabric coördinates

[vel, pvecs] = christof(n, cc, den); % in fabric coördinates

pvecs = rotall'*pvecs; % in global directions
fvec = pvecs(:,2);
fvec = fvec(2:3)./norm(fvec(2:3));
%% Outputs
fastaz = r2d(atan2(fvec(1),fvec(2)));
DinvV = 1./vel(1) - 1./vel(2);
Vav = 0.5*(vel(1) + vel(2));
anis = 100*(vel(2)-vel(1))/(2*Vav);
end