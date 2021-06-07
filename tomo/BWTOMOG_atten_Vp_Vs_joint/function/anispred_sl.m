function [phipred,dTpred,SIpred] = anispred_sl(phitrue,dTtrue,period,betas,plotopt)
% [phipred,dTpred,SIpred] = anispred_sl(phitrue,dTtrue,period,betas,plotopt)
%
% Adapted from Silver and Long 2011, this script analytically calculates
% the theoretical anisotropy parameters experienced over N layers by
% incoming rays with polarisations "betas". 
% N.B. ray parameter is ignored. 
% 
% INPUTS:
% phitrue = [1 x N] vector of layer fast azimuths ordered from deepest
%           layer to shallowest (i.e. deepest first)
% dTtrue  = [1 x N] vector of layer delay times ordered from deepest
%           layer to shallowest (i.e. deepest first)
% period  = period of the seismic waves (seconds)
% betas   = vector of arrival azimuths (i.e. SKS polarisations) in degrees
%           (default is to do 0:1:359 degrees if no argument given)
%
% OUTPUTS
% phipred = predicted net values of phi for each arrival in betas
% dTpred = predicted net values of dT for each arrival in betas
% SIpred = predicted net values of splitting intensity for each arrival in betas
%
% Written by Zach Eilon, 2013


if nargin<4
    betas=0:1:359;
end
if nargin <5
    plotopt=0;
end

omega = 2*pi/period;

nlayers = length(phitrue);
if length(dTtrue)~=nlayers, error('phitrue and dTtrue must have same number of elements'); end
nbetas = length(betas); 


ars = zeros(nbetas,1); % alpha_radial
ats = zeros(nbetas,1); % alpha_transverse
Cc = zeros(nbetas,1);
Cs = zeros(nbetas,1);

phipred = zeros(nbetas,1); % predicted net fast azimuth
dTpred  = zeros(nbetas,1); % predicted net time delay
SIpred  = zeros(nbetas,1); % predicted net splitting intensity

for ip=1:nbetas
    as = d2r(2*(phitrue - betas(ip))); % alphas
    as = mod(as+pi,2*pi)-pi;
    ths = omega.*dTtrue/2; % thetas

    S = 1; % S = product over layers of cos(theta)
    for ii=1:nlayers;        S = S*cos(ths(ii));    end
    
    arsum = 0;
    for ii = 1:nlayers-1
        for jj = ii+1:nlayers
           arsum = arsum + tan(ths(ii))*tan(ths(jj))*cos(as(ii) - as(jj));
        end
    end
	atsum = 0;
    for ii = 1:nlayers-1
        for jj = ii+1:nlayers
           atsum = atsum + tan(ths(ii))*tan(ths(jj))*sin(as(ii) - as(jj));
        end
    end
    Ccsum = 0;
    for ii = 1:nlayers
        Ccsum = Ccsum + tan(ths(ii))*cos(as(ii));
    end
    Cssum = 0;
    for ii = 1:nlayers
        Cssum = Cssum + tan(ths(ii))*sin(as(ii));
    end
    
    ars(ip) = S*(1-arsum);
    ats(ip) = S*atsum;
    Cc(ip)  = S*Ccsum;
    Cs(ip)  = S*Cssum;
end

phipred = zeros(nbetas,1);
dTpred  = zeros(nbetas,1);
for ip=1:nbetas
    thetap = atan(sqrt((ars(ip)*ats(ip) + Cs(ip)*Cc(ip))^2 + (ats(ip)^2 + Cs(ip)^2)^2)/(ars(ip)*Cs(ip) - ats(ip)*Cs(ip)));
    sinap = ((ats(ip)^2 + Cs(ip)^2)/(ars(ip)*Cs(ip) - (ats(ip)^2)*Cc(ip)))/tan(thetap);
    cosap = ((ars(ip)*ats(ip) + Cs(ip)*Cc(ip))/(ars(ip)*Cs(ip) - (ats(ip)^2)*Cc(ip)))/tan(thetap);
    
%     sinap = ats(ip)^2 + Cs(ip)^2;
%     cosap = ars(ip)*ats(ip) + Cs(ip)*Cc(ip);
    
    alphap = atan2(sinap/tan(thetap),cosap/tan(thetap));
    phipred(ip) = betas(ip) + 0.5*mod(r2d(alphap),360);
    phipred(ip) = mod(phipred(ip)+90,180)-90;
    dTpred(ip) = 2 * abs(thetap) /omega;
    SIpred(ip) = (ats(ip)^2 + Cs(ip)*Cc(ip))/(ars(ip)*Cs(ip) - (ats(ip)^2)*Cc(ip));
end
    
% %%%
if plotopt==1
figure(1)
subplot(211)
plot(betas,phipred,'-')
set(gca,'XLim',[0 360],'YLim',[-90,90],'YTick',[-80:20:80],'YGrid','on');
subplot(212)
plot(dTpred,'-')
ylim([0 4])
end

disp('NB problem with the SIprid part of this. Gives SI 1/3 of expected')

end
    
    
    