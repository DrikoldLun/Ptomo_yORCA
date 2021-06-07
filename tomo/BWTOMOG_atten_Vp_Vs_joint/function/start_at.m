function [ start_index ] = start_at( xx,min_x )
%[ start_index ] = start_at( xx,min_x )
%   function to find the index to begin iteration partway through a list of
%   non-consecutive indices. E.g. You are halfway through a list of events:
%   [1 2 5 8 15], and you know you've done up to orid 4, this function will
%   return the index 3.
% 
% the index will start you in the list xx at the first element that is equal
% to or greater then min_x

start_index = ceil(zerof(xx,1:length(xx),min_x));

if min_x < min(xx), start_index=1; end
if min_x > max(xx), error('Start value exceeds maximum in series'); end

end

