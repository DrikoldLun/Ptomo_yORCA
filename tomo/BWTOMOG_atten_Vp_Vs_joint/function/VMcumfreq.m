function [FKth] = VMcumfreq(th,K,mu)
%function [FKth] = VMcumfreq(th,K,I0K)
%Approximation to true interal for von Mises CDF
% th=0.1;
I0K=besseli(0,K);
dph=0.00001;
x=[-pi:dph:th];
FKth = trapz(x,exp(K*cos(x)));
FKth= FKth/(2*pi*I0K);


% sumth=0;
% for j=1:1000
%     sumth = sumth + besseli(j,K)*sin(j*th)/j;
% end
% sumth0=0;
% for j=1:1000
%     sumth0 = sumth0 + besseli(j,K)*sin(j*th0)/j;
% end
% FKth=(1/2*pi)*((th - th0) + (2/I0K)*(sumth-sumth0))











end

