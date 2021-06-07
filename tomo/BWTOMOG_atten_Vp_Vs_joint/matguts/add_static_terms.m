function [ d,G ] = add_static_terms( d,G,data,model )
% subtract event and station terms from d and G matrix
% positive e_stat or s_stat will decrease predicted dt/tstar.
% therefore these are additional delays - positive means slow structure
% beneath just that station.
% NB adds zero rows for orids or stas are not in the data
%
% model parameter vector will be [mdq estatic sstatic]

%% event terms
rayorids = data.ray.orid;
orids = data.evt.orid;
nevts = data.evt.nevts; 

e_d = zeros(size(d)); % event term for each datum (ray/observation)
e_G = spalloc(size(G,1),nevts,nevts);
for ie = 1:nevts
    ind = rayorids==orids(ie);
    e_d(ind) = model.estatic(ie);
    e_G(:,ie) = -ind;
end

%% station terms
raystas = data.ray.sta_num;
stas = data.stn.num;
nstas = data.stn.nstas;
phs = data.ray.ph;

% find out if p and s or just one
uphs = unique(phs);
nph = length(uphs);


s_d = zeros(size(d)); % station term for each datum (ray/observation)
s_G = spalloc(size(G,1),nstas*nph,nstas*nph);
for ip = 1:nph
    addsta = (ip-1)*nstas;
    for is = 1:nstas
        ind = raystas==stas(is) & phs==uphs(ip);
        s_d(ind) = model.sstatic(is + addsta);
        s_G(:,is + addsta) = -ind;
    end
end

d = d - e_d - s_d; % positive e_stat or s_stat will decrease t_pred
G = [G,e_G,s_G];

end
