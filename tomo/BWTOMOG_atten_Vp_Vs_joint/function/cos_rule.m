function [ a ] = cos_rule( b,c,A )
%[ a ] = cos_rule( b,c,A )
%   applies the cosine rule to get the third side of the triangle (a) from
%   the other two lengths (b,c) and the angle between them (A) in degrees
% a^2 = b^2 + c^2 - 2*b*c*cos(A)

a = sqrt( b.^2 + c.^2 - 2*b.*c.*cosd(A) );



end

