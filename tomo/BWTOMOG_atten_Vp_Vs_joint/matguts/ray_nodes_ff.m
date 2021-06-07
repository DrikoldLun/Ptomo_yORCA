function [n_indx, n_vals ] = ray_nodes_ff(rayxyzdr,r_cdr,rmax,par,p,cf,trdist)
% [n_indx, n_vals ] = ray_nodes_ff(rayxyzdr,r_cdr,rmax,par,p,cf,trdist)
% 

% first pass - find nodes nearby ray
kdr = rayxyzdr(2,4)-rayxyzdr(1,4);
rmax = sqrt( (4*kdr)^2 + rmax^2);
nn = 2:4:length(rayxyzdr(:,1))-1;
rayxyz2 = [ rayxyzdr(nn,1) rayxyzdr(nn,2) rayxyzdr(nn,3) ];
 
nn=1;
n_indx = zeros(9000,1);
for qq = 1:length(rayxyz2(:,1))
    dists = ( (par.mx - rayxyz2(qq,1)).^2 ...
            + (par.my - rayxyz2(qq,2)).^2 ...
            + (par.mz - rayxyz2(qq,3)).^2 ).^(1/2);
    ind1 = find(dists<rmax);
    n_indx(nn:(nn+length(ind1)-1)) = ind1;
    nn = nn+length(ind1);
end

n_indx = n_indx(n_indx>0);
n_indx = unique(n_indx);

l_nindx = length(n_indx);

% calc Rn for all nearby nodes
r_scale = (6371-rayxyzdr(:,3))/6371;
n_rn = zeros(l_nindx,1);
n_dr = zeros(l_nindx,1);
n_cdr = zeros(l_nindx,1);
for qq = 1:length(n_indx)
   dists = ( (r_scale.*(rayxyzdr(:,1)-par.mx(n_indx(qq)))).^2 ...
           + (r_scale.*(rayxyzdr(:,2)-par.my(n_indx(qq)))).^2 ...
           + (rayxyzdr(:,3)-par.mz(n_indx(qq))).^2 ).^(1/2);   
            ind1 = find(dists==min(dists));
            n_rn(qq) = dists(ind1);
            % these last two are actually equal... should they be?
            n_dr(qq) = rayxyzdr(ind1,4);
            n_cdr(qq) = r_cdr(ind1);
end

% calc mrdist
mrdist = max(n_cdr)-min(n_cdr);

% function to calc sensitivity vals
[ n_indx, n_vals ] = f1_dtdm(n_indx,n_rn,n_dr,n_cdr,mrdist,par,p,cf,trdist);

end
