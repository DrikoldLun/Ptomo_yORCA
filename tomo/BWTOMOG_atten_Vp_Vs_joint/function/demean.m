function [ dV,Vav ] = demean( V )
%[ dV , Vav] = demean( V )
% quick function to remove mean from columns of V

Vav = mean(V);
dV = detrend(V,'constant');


end

