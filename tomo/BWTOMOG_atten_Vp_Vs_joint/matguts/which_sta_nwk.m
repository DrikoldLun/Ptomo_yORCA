function nwk = which_sta_nwk(stainfo,sta,slat,slon)
% nwk = which_sta_nwk(stainfo,sta,slat,slon)
%
% take stainfo structure, and for a certain named station, with lat and
% long, find the network

isname = find(strcmp(stainfo.stas,sta));

if length(isname)>1
    % find closest station
    isname = isname(mindex(distance(slat,slon,stainfo.slats(isname),stainfo.slons(isname))));
end

nwk = stainfo.nwk{isname};

end

