function [minang,maxang] = angrange(depth,phase,vSsurf)
% [minang,maxang] = angrange(depth,phase,vSsurf)
% function to give the
% angular window for which a given phase arrives at surface station inputs
% are: depth of event (in km), phase required (string with quotes) and the
% surface velocity (optional)
%
% Written by Zach Eilon, 2012


% depth=20;
% phase='SKS';

% Get near-surface velocity from input, or use default value
if nargin==3
    vSsurf=vSsurf;
else
    vSsurf=3; %default value
end

%test over all angular distances
testdeg=zeros(180,1);
for i=1:length(testdeg)
    tt=tauptime([],depth,phase,'deg',i);
    if isempty(tt)==0
        testdeg(i)=tt.time;
    end
end
minang=find(testdeg>0,1,'first');
maxang=find(testdeg>0,1,'last');

% Test income angles for Scrtical criterion
inc_min=90;
while inc_min > 32
tt=tauppath([],depth,phase,'deg',minang);
inc_min=asind(vSsurf*tt.rayParam/6371); % according to ? = R*sin(i)/ Vs
if inc_min > 32
    minang=minang+1;
end
end
fprintf('Minang is %3i, with an incidence angle at surface of %4.1f (degrees) \n',minang,inc_min);

inc_max=90;
while inc_max > 32
tt=tauppath([],depth,phase,'deg',maxang);
inc_max=asind(vSsurf*tt.rayParam/6371); % according to ? = R*sin(i)/ Vs
if inc_max > 32
    maxang=maxang-1;
end
end
fprintf('Maxang is %3i, with an incidence angle at surface of %4.1f (degrees) \n',maxang,inc_max);


    
    
    
    