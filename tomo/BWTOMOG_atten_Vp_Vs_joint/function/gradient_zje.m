function [ dFdX ] = gradient_zje( F,X )
%[ FX ] = gradient_zje( F,X )
%   function to calculate gradient of a series with uneven spacing in the
%   dimension along which you are differentating.
% NB RELIES UPON CLOSE SPACING - will fit cubic spline, but still, gradient
% only calculated based on (f(x)-f(x+dx))/dx
N = length(F);
if length(X)~=N; error('F and X must be same length'); end
df_dx = diff(F)./diff(X);
midx = 0.5*(X(1:N-1) + X(2:N));
dFdX = interp1(midx,df_dx,X,'pchip','extrap');
end

