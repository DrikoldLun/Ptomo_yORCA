function [ n_vals ] = f1_dtdm_wt4(dists,n_vel,n_vals1,n_cdr)

n_vals = zeros(size(n_vals1));
for idist = 2:length(dists)-1
    % work out indices of nodes at this distance and on either side
    ind1 = find(n_cdr==dists(idist));
    ind2 = find(n_cdr==dists(idist-1));
    indx = union(ind1,ind2);
    ind3 = find(n_cdr==dists(idist+1));
    indx = union(indx,ind3);
    
    rl = dists(idist+1)-dists(idist-1);
    % ZE: this must be the normalisation??
    total = rl/mean(n_vel(indx)); % time
    
    n_vals(indx) = n_vals1(indx)/sum(n_vals1(indx)); % normalised vals
    n_vals(indx) = -total*n_vals(indx); % this nvals is then in units of seconds*K
                              % K is then unitless, and mparms are fractional slowness...
end

end

