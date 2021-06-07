function [ N ] = unhandshake( npair )
%[ N ] = unhandshake( npair )
% from npair handshakes between "N" persons, i.e. npair = N(N-1)/2
% back out the number of persons (simple quadratic formula)

N = (1+sqrt(1 + 8*npair))/2;

end

