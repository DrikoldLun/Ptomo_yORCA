function [gcarc] = angsep(lat1,long1,lat2,long2)
% [gcarc] = angsep(lat1,long1,lat2,long2)
% This function calculates only the angular distance between two points on
% the Earth's surface, given geographical coordinates.
% It uses the "ellipsep" function, but only outputs the gcarc
%
% Written by Zach Eilon, 2013



if length(lat2)==length(long2) && length(lat1)==1 && length(long1)==1 
gcarc=zeros(size(lat2));   
for i=1:length(lat2)
[foraz,backaz,gcdist,gcarci] = ellipsep (lat1,long1,lat2(i),long2(i));
gcarc(i)=gcarci;
end
else
    error('angsep','Size of input arguments is not consistent')
end