function [ n_indx, n_vals ] = f1_dtdm(n_indx,n_rn,n_dr,n_cdr,mrdist,par,p,cf,trdist)

n_vel = par.mvav(n_indx);

wl = n_vel./cf;  % wavelength (=velocity/freq)

ndh = 0.5*(par.mdx(n_indx) + par.mdy(n_indx)) ;
dh_rn = ndh.*( ((6371-par.mz(n_indx))/6371) );

RF1 =  ((wl.*(n_dr.*(trdist-n_dr)))/trdist).^(1/2); % Fresnel zone radius
RF1_min = 0.9*dh_rn;
too_small = RF1./RF1_min;
RF1(too_small<1) = RF1_min(too_small<1);

RF1_frac = n_rn./RF1;
indx = find(RF1_frac<1.3);
n_indx = n_indx(indx);
n_rn = n_rn(indx);
RF1 = RF1(indx);
n_vel = n_vel(indx);
n_dr = n_dr(indx);
n_cdr = n_cdr(indx);
dh_rn = dh_rn(indx);

n_dr2 = 0.33*(n_dr.^(1/2));
wt = zeros(length(n_indx),11);
wt(:,1) = (1-((n_rn./(1.4*RF1)).^2)).*sin(pi*((n_rn+n_dr2)./(n_dr2+RF1)).^2);
wt(:,2) = (1-(((n_rn+0.07*dh_rn)./(1.4*RF1)).^2)).*sin(pi*((n_rn+0.07*dh_rn+n_dr2)./(n_dr2+RF1)).^2);
wt(:,3) = (1-(((n_rn-0.07*dh_rn)./(1.4*RF1)).^2)).*sin(pi*((n_rn-0.07*dh_rn+n_dr2)./(n_dr2+RF1)).^2);
wt(:,4) = (1-(((n_rn+0.15*dh_rn)./(1.4*RF1)).^2)).*sin(pi*((n_rn+0.15*dh_rn+n_dr2)./(n_dr2+RF1)).^2);
wt(:,5) = (1-(((n_rn-0.15*dh_rn)./(1.4*RF1)).^2)).*sin(pi*((n_rn-0.15*dh_rn+n_dr2)./(n_dr2+RF1)).^2);
wt(:,6) = (1-(((n_rn+0.24*dh_rn)./(1.4*RF1)).^2)).*sin(pi*((n_rn+0.24*dh_rn+n_dr2)./(n_dr2+RF1)).^2);
wt(:,7) = (1-(((n_rn-0.24*dh_rn)./(1.4*RF1)).^2)).*sin(pi*((n_rn-0.24*dh_rn+n_dr2)./(n_dr2+RF1)).^2);
wt(:,8) = (1-(((n_rn+0.33*dh_rn)./(1.4*RF1)).^2)).*sin(pi*((n_rn+0.33*dh_rn+n_dr2)./(n_dr2+RF1)).^2);
wt(:,9) = (1-(((n_rn-0.33*dh_rn)./(1.4*RF1)).^2)).*sin(pi*((n_rn-0.33*dh_rn+n_dr2)./(n_dr2+RF1)).^2);
wt(:,10) = (1-(((n_rn+0.44*dh_rn)./(1.4*RF1)).^2)).*sin(pi*((n_rn+0.44*dh_rn+n_dr2)./(n_dr2+RF1)).^2);
wt(:,11) = (1-(((n_rn-0.44*dh_rn)./(1.4*RF1)).^2)).*sin(pi*((n_rn-0.44*dh_rn+n_dr2)./(n_dr2+RF1)).^2);
wt = real(wt);
wt(wt<0) = 0;
wt = mean(wt,2);

n_vals = wt;

dists = unique(sort(n_cdr));
n_vals = f1_dtdm_wt4_zemod(dists,n_vel,n_vals,n_cdr);

n_indx = n_indx(abs(n_vals)>0);
n_vals = n_vals(abs(n_vals)>0);

end

