function [ data ] = crustal_correction( data,par )
% [ data ] = crustal_correction( data,par )
%   performs crustal correction at each station
%   at this point, assumes a uniform crustal S-wave velocity
% 
% Correcting for rays spending RELATIVELY more/less time in the crust
% compared to arrivals at other stations from the same earthquakes
% a POSITIVE time corresponds to more time in the crust (thicker than av)
% a NEGATIVE time corresponds to less time in crust (thinner than av)
%  thus, these corrections should be SUBTRACTED from measured dT to give
%  the contribution to the dT that arises from the mantle. 



% %% Check crust_corr.mat has been made, and give option to re-do
% if exist('crust_corr.mat','file')==2
%     if par.redo_crust==1
%         redo = input('Re-do the crustal correction calc? y/n [n] ','s');
%     else
%         redo = 'n';
%     end
% else
%     redo = 'y';
% end
% if strcmp(redo,'y')
%     run('calc_crustcor_plot.m');
% end
% if exist('crust_corr.mat','file')~=2
%     fprintf('crust_corr still not there, something is wrong... stopping\n')
%     return
% end
% cc = load('crust_corr.mat');
    

if par.crust_corr % check to see if we are doing crustal correction!
   
ray = data.ray;     
stn = data.stn;
% [~,icc,~] = intersect(cc.stas,stn.sta);
ray.ccorr = zeros( size(ray.p) );
unique_orid = unique( ray.orid );   
norid = length( unique_orid );   
% set up differential correction time.
% this corrects for a ray spending RELATIVELY more/less time in the crust
% compared to arrivals at other stations from the same earthquakes
% a POSITIVE time corresponds to more time in the crust (thicker than av)
% a NEGATIVE time corresponds to less time in crust (thinner than av)
%  thus, these corrections should be SUBTRACTED from measured dT to give
%  the contribution to the dT that arises from the mantle. 
dtcor = zeros(size(ray.d)); 

%% NO correction if doing Q inversion - no prior for what crustal Q should be like
if par.t_ts == 2
    data.ray = ray;
    warning('No crustal correction for Q inversion')
    return
end

% average properties of crust - can change these but probably pretty good
vpav = 6.56;
vpvsav = 1.8;
vum = [8 4.4]; % p, s
stn.vc(:,1) = vpav*ones(stn.nstas,1);
stn.vpvs(isnan(stn.vpvs)) = vpvsav; % use average value if vpvs not known
stn.vc(:,2) = vpav./stn.vpvs;


if par.PS == 1
    vind = 1;
elseif par.PS == 2 
    vind = 2;
end


fprintf('>  Computing crust corrections... \n');
for ie = 1:norid                	
	
    indx = find( ray.orid == unique_orid(ie) );
    sta_num = ray.sta_num( indx );  
    if isfield(ray,'ph'), vind = unique(ray.ph(indx)); end
    
    % elevation time correction (assume pure vertical)
    Telv = stn.elv(sta_num)./stn.vc(sta_num,vind);

    % moho thickness time correction
    % moho depth from sea level...
    mohz = stn.moh(sta_num) - stn.elv(sta_num); m_mohz = nanmean(mohz);
    % set up time correction
    dTmoh = zeros(length(indx),1);
    for is = 1:length(indx)
        snum = sta_num(is);
        if mohz(is)==0 || isnan(mohz(is))
            dTmoh(is) = nan;
        else
            % only vertical travel time
%             dTmoh(is) = (mohz(is)-m_mohz)*(1./stn.vc(snum,vind) - 1./vum(vind));
            % including inc... verttime / cosine(inc_av) - use average incidence angle in crust 
            cinc = rayp2inc(ray.pd(indx(is)),stn.vc(snum,vind),6371-mohz(is)/2); % use half crustal thickness as depth
            minc = rayp2inc(ray.pd(indx(is)),vum(vind),        6371-mohz(is)); % at moho
            dTmoh(is) = (mohz(is)-m_mohz)*(1./stn.vc(snum,vind)./cosd(cinc) - 1./vum(vind)./cosd(minc));
        end
    end
    
    % make sure de-meaned
    dTmoh = dTmoh - nanmean(dTmoh);
    % set correction for stations without known moho to 0
    dTmoh(isnan(dTmoh)) = 0;
    dTelv = Telv - mean(Telv);
    
    ray.ccorr(indx) = dTmoh + dTelv;
end

data.ray = ray;


else 
    fprintf('>  NO crust corrections \n');
    data.ray.ccorr = zeros(data.ray.nrays,1);
end

end

