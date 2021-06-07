function [LATOUT, LONOUT] = reckon_km(LAT, LON, ARCLEN, AZ)
% [LATOUT, LONOUT] = reckon_km(LAT, LON, ARCLEN, AZ)
%   This is just the reckon function but with the Earth's ellipsoid input
%   by defult, such that the arclen input can be entered in km.

[LATOUT,LONOUT] = reckon(LAT,LON,ARCLEN,AZ,[6378.137, 0.08181919]);


end

