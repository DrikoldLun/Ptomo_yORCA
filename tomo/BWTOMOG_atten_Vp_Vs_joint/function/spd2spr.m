function [ p_spr ] = spd2spr( p_spd )
%[ p_spr ] = spd2spr( p_spd )
%   Convert a slowness in seconds per degree to seconds per radian
%  uses Earth radius of 6371
%  p_spr = (180/pi)*p_spd/Re

p_spr = r2d(p_spd)/6371;


end

