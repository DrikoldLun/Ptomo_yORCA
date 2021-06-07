function plot_crust_corr(data,par,saveres)
% plot_data(data,par,saveres)

if nargin < 3
    saveres = false;
end

figure(32), clf, hold on
% mkfig_CascMAP

lonlims = par.plot_lonlims;
latlims = par.plot_latlims;
plotsize = 1200;
hlim = [20 65];
dlim = [-2 2];


scale = 0.03; % distance multiplier


indx = abs(data.ray.ccorr) < 4; 
dTcrust = data.ray.ccorr(indx);
baz = data.ray.baz(indx);
slat = data.ray.stalat(indx);
slon = data.ray.stalon(indx);

% crustal velocities etc.
vpav = 6.56;
vpvsav = 1.8;
vsav = vpav./vpvsav;


if par.PS == 1
    psstr = 'P';
    if  ~isfield(data.ray,'ph'), data.ray.ph = ones(data.ray.nrays,1); end
elseif par.PS == 2
    psstr = 'S';
    if  ~isfield(data.ray,'ph'), data.ray.ph = 2*ones(data.ray.nrays,1); end
end

% compute incidence angles
rayps = [1:1:15]';
incs(:,1) = rayp2inc(rayps,vpav,6371);
incs(:,2) = rayp2inc(rayps,vsav,6371);

rincs = nan(data.ray.nrays,1);
rincs(data.ray.ph==1) = interp1(rayps,incs(:,1),data.ray.pd(data.ray.ph==1));
rincs(data.ray.ph==2) = interp1(rayps,incs(:,2),data.ray.pd(data.ray.ph==2));

% compute horiz distances
xdists = data.stn.moh(data.ray.sta_num).*tand(rincs);

% compute pierce points
plo = data.stn.lon(data.ray.sta_num) + scale*xdists.*sind(baz);
pla = data.stn.lat(data.ray.sta_num) + scale*xdists.*cosd(baz);

    

%% plot uncorrected times on Map
set(gcf,'position',[100 200 plotsize/plot_size_ratio(lonlims,latlims) plotsize])
set(gca,'xlim',lonlims,'ylim',latlims)
hold on

%% plot the data
% plot station mohos
for is = 1:data.stn.nstas
    if isnan(data.stn.moh(is)), continue; end
    plot(data.stn.lon(is),data.stn.lat(is),'o','Markersize',10,...
        'MarkerFacecolor',colour_get(data.stn.moh(is),hlim(2),hlim(1),jet),...
        'MarkerEdgeColor',0.6*[1 1 1])
end

% plot P arrivals
isp = data.ray.ph==1;
scatter(plo(isp),pla(isp),10,dTcrust(isp),'v','filled')
% plot S arrivals
iss = data.ray.ph==2;
scatter(plo(iss),pla(iss),10,dTcrust(iss),'o','filled')


hcb = colorbar('peer',gca);
colormap(parula)
caxis(dlim)
set(get(hcb,'Ylabel'),'string',['Crustal correction (s)'],...
    'Fontsize',24,'FontWeight','bold','interpreter','Latex');
set(get(hcb,'Ylabel'),'position',get(get(hcb,'Ylabel'),'position')+[0.4 0 0])


% % plot on a key of weights
% wt_eg = [1, 5, 10];
% sd_eg = (10*wt_eg).^-0.5;
% wd_eg= 0.5./sd_eg;
% lon_key = -131;
% lat_key = 41;
% hold on
% for ii = 1:length(wt_eg)
%     plot(lon_key + [0, scale],(lat_key + 0.1*ii)*[1 1],'Linewidth',wd_eg(ii),'color','k')
%     text(lon_key + 1.5*scale,   lat_key + 0.1*ii,num2str(wt_eg(ii)),'Fontsize',15,'FontWeight','bold')
% end
% plot key with datatype
% hold on
% m_text(lon_key,-10.85,['\textbf{',secstrs{ic},'}'],'Fontsize',22,'FontWeight','bold','interpreter','Latex')

% save
if saveres
    ofile = ['map_data_',datstr,'_',psstr];
    save2pdf(32,ofile,par.figdir);
end


end

