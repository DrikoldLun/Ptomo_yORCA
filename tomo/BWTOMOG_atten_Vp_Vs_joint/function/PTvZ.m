function [  P,T,rho ] = PTvZ( Z,moho )
%[ P,T ] = PTvZ( Z,moho )
%   Function to calculate expected temperature and pressure at a given
%   depth within the Earth. Only works down to transition zone
% INPUTS
%  Z    = a vector of depths (km)
%  moho = moho depth (km)
% OUTPUTS
%  P    = a vector of pressures (GPa)
%  T    = a vector of temperatures (K)
%  rho  = a vector of densities (g/cc)


P = zeros(size(Z));
T = zeros(size(Z));
rho = zeros(size(Z));

zw  = 0.000;
zc1 = 13;
zmo = moho;
zli = 120;

g = 9.81;


frho  = [1.02  1.02 2.6  2.6  2.9  2.9  3.2  3.1 3.5];
zrho =  [0     zw   zw   zc1  zc1  zmo  zmo  zli 440];

fT = [0 30*(zmo) (0.6*(440-zmo) + 30*zmo)]+283;
zT = [0 zmo      440                     ];

for iz = 1:length(Z)
    rho(iz) = linterp(zrho,frho,Z(iz));
    kk = sum(zrho < Z(iz));
    if kk<2, P(iz)=0; T(iz)=fT(1); continue, end
    for ik = 1:kk-1
        ztop = zrho(ik);
        zbot = zrho(ik+1);
        if zbot==ztop, continue; end
        mean_rho = 0.5*(linterp(zrho,frho,zbot-0.001)+linterp(zrho,frho,ztop+0.001));
        P(iz) = P(iz) + g*mean_rho*(zbot-ztop)*10^6;
    end
    P(iz) = P(iz) + g*(Z(iz)-zbot)*0.5*(linterp(zrho,frho,Z(iz)-0.001)+linterp(zrho,frho,zbot+0.001))*10^6;
    T(iz) = linterp(zT,fT,Z(iz));
end


P = P/1e9;

end