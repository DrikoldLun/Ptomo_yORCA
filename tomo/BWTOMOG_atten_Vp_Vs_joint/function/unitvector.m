function [ v_norm ] = unitvector( v )
% [ v_norm ] = unitvector( v )
%   function to make vector(s) v of unit length

v_norm = v./norm(v);


end

