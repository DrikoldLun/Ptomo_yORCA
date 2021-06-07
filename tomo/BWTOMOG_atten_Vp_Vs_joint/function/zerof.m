function [ xx0 ] = zerof( yy,xx,level )
% [ xx0 ] = zerof( yy,xx,level )
% For a function yy of xx, find the xx values correspoinding to all the
% (linear interpolated) crossing points where yy goes to zero


N = length(yy);

if nargin < 2
    xx=1:N;
end
if nargin==3
    yy = yy-level; % shift to zero at level
end
    
    
if length(xx)~=N, error('xx and yy must be same length'); end
xx0=zeros(N,1);

for ii = 1:N-1
    if yy(ii)==0 % check if identically zero
        xx0(ii)=xx(ii);
    elseif sign(yy(ii))~=sign(yy(ii+1))
        if yy(ii+1)==0, continue; end % if next one is identically zero, wait till next time
        % interpolate: x = x1 + dx(a/(a+b))
    xx0(ii)= xx(ii) + (xx(ii+1)-xx(ii))*(abs(yy(ii))/(abs(yy(ii))+abs(yy(ii+1))));
    end
end
if yy(end)==0
    xx0(end)=N;
end

xx0 = xx0(xx0~=0);

end

