function [ xx ] = ceil_level( x, y )
%[xx] = ceil_level(x,level)
% ceils the number/vector x to the nearest y (the level)
%
% Written by Zach Eilon, 2012

xx = ceil(x./y).*y;
end

