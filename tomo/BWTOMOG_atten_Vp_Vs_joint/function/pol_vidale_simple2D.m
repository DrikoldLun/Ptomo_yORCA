function [ polest,polstrength ] = pol_vidale_simple2D( datN,datE )
%[ polest ] = pol_vidale_simple2D( datN,datE )
%   Simple function to estimate the polarisation of a signal on the N and E
%   components
% 
% polest is in degrees

%% Measure initial polarisation following Vidale 1986
u = datE + 1i*hilbert(datE);
v = datN + 1i*hilbert(datN);
CC = [u'*conj(u) u'*conj(v)
      v'*conj(u) v'*conj(v)];
% C = [datN'*datN, datN'*datE; datE'*datN, datE'*datE];
[V,D] = eigs(CC);
V = real(V);
D = real(diag(D));

S = V(1,D==max(D));
C = V(2,D==max(D));
if S>=0 && C>0
    polest=atan(S/C); % in radians
elseif S<0 && C>0
    polest=atan(S/C) + 2*pi; % in radians
elseif C<=0
    polest=atan(S/C)+pi; % in radians
end

polest = r2d(polest);
polstrength = 1 - min(D)/max(D);

end

