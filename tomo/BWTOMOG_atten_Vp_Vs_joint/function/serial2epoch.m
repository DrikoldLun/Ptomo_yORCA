function [epochtime] = serial2epoch(serialtime)
%[epochtime] = serial2epoch(serialtime)
%   Convert serial time (days since 01-01-0000) to epochal time (seconds
%   since 01-01-1970

epoch0 = datenum(1970,1,1);

epochtime = (serialtime-epoch0)*24*60*60;

end

