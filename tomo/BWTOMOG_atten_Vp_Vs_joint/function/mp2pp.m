function [xx] = mp2pp( x )
%function [xx] = mp2pp( x ) 
%converts a vector with elements (in radians) between 0 and 2pi to the range
%-pi to +pi
x=mod(x,2*pi);
x(find(x > pi)) = x(find(x > pi)) - 2*pi;
xx=x;
end

