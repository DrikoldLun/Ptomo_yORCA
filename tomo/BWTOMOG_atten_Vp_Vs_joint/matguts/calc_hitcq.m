function [ par ] = calc_hitcq( data,par,K )
% [ par ] = hitmap( data,par,K )
%  function to calculate how many rays go through each voxel in the model

% maxNrays = 5;
maxNrays = 10;
sectors  = 6;
if par.PS == 1
    Nph = 1; phm = 0;
elseif par.PS == 2
    Nph = 1; phm = -1; % account for the fact that now .ph=2 but one field
elseif par.PS == 3
    Nph = 2; phm = 0;
end

    

hitc = zeros(par.nmodel,sectors,Nph);
for ii = 1:length(K.n_indx)
    % find hexant (1 is 0-60, 2 is 60-120 etc.)
    sect = ceil(mod(data.ray.baz(ii),360)/(360/sectors));
    if sect==0, sect=1; end % in pathological case... 
    % assign to hexant
    hitc(K.n_indx{ii},sect,data.ray.ph(ii)+phm) = hitc(K.n_indx{ii},sect,data.ray.ph(ii)+phm) + 1; % simple count of passing box
%     hitc(K.n_indx{ii},sect,data.ray.ph(ii)) ...
%         = hitc(K.n_indx{ii},sect,data.ray.ph(ii)) + K.n_vals{ii}./sum(K.n_vals{ii}); % account for sensitivity
end


temp = hitc;
temp(hitc>maxNrays)=maxNrays; % cap at minNrays
hitq = squeeze(sum(temp,2))./(maxNrays*sectors); % quality = hits out of a possible maximum of cap*sectors

par.hitc=hitc;
par.hitq=hitq;

end

