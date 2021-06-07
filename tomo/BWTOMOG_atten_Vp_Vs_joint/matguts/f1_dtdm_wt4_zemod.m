function [ n_vals ] = f1_dtdm_wt4_zemod(dists,n_vel,n_vals1,n_cdr)

n_vals = zeros(size(n_vals1));
% tottt = 0;
for idist = 2:length(dists)%-1
    % work out indices of nodes at this distance and on either side
    ind1 = find(n_cdr==dists(idist));
    ind2 = find(n_cdr==dists(idist-1));
    indx = union(ind1,ind2);
% % ZE comment - changed so we only consider segments between each pair of
% node layers, not three layers at a time...
%     ind3 = find(n_cdr==dists(idist+1)); 
%     indx = union(indx,ind3); 
    
    rl = dists(idist)-dists(idist-1);
    % ZE: time normalisation - so ray integrates to correct value
    total = rl/mean(n_vel(indx)); % time
%     tottt = tottt+total;
    
% % ZE comment - calc. this separately, as was overwritng n_vals(indx)
% things before.
    n_vals_norm = n_vals1(indx)/sum(n_vals1(indx)); % normalised vals
    n_vals(indx) = n_vals(indx) -total*n_vals_norm; % this nvals is then in units of seconds*K
                       % K is then unitless, and mparms are fractional slowness...
end
% tottt
% sum(n_vals)

end

