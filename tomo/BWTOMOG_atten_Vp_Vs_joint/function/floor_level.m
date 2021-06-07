function [ xx ] = floor_level( x, y )
%[xx] = floor_level(x,level)
% floors the number/vector x to the nearest y (the level)
xx = floor(x./y).*y;
end

