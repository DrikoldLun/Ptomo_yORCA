function [serialtime] = epoch2serial(epochtime)
%[serialtime] = epoch2serial(epochtime)
%   Convert epochal time (seconds since 01-01-1970 to serial time (days
%   since 01-01-0000)

epoch0 = datenum(1970,1,1);

serialtime = epoch0 + epochtime./(24*60*60);

end

