function [ theta_ ] = m90_p90( theta )
%[ theta_ ] = m90_p90( theta )
%   function to take any angle and put it onto the range -90 < theta <= +90
% e.g. -130 ==>  50
% e.g.   91 ==> -89
% e.g.  207 ==>  27

theta_ = mod(theta+90,180) - 90;
if theta_ == -90, theta_=90; end

end

