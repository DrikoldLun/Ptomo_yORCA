function plot_data(data,par,saveres)
% plot_data(data,par,saveres)

figure(23), clf, hold on

lonlims = par.plot_lonlims;
latlims = par.plot_latlims;
plotsize = 800;
dlim = [-0.5 0.5];


scale = 0.4; % length of lines


indx = abs(data.ray.d) < 4; 
d = data.ray.d(indx);
baz = data.ray.baz(indx);
sd = data.ray.sd(indx);
slat = data.ray.stalat(indx);
slon = data.ray.stalon(indx);

if par.PS == 1
    psstr = 'P';
elseif par.PS == 2
    psstr = 'S';
elseif par.PS == 3
    psstr = 'P+S';
end
if par.t_ts == 1;
    datstr = 'dT';
    lbstr = '\delta T';
    cmap = jet;
    dlim = 1.25*dlim;
elseif par.t_ts == 2;
    datstr = 'dtstar';
    lbstr = '\Delta t^*';
    cmap = parula;
end

 
    

%% plot uncorrected times on Map
set(gcf,'position',[200 300 plotsize/plot_size_ratio(lonlims,latlims) plotsize])
set(gca,'xlim',lonlims,'ylim',latlims)


for ia=1:sum(indx)

% plot the data
p1la = slat(ia);
p1lo = slon(ia);
[p2la,p2lo] = reckon(p1la,p1lo,scale,baz(ia));
% plot line toward baz, coloured by tstar
plot([p1lo;p2lo],[p1la;p2la],'LineWidth',0.2/sd(ia),...
    'color',colour_get(d(ia),dlim(2),dlim(1),cmap))


end % loop on arrivals

hcb = colorbar('peer',gca);
colormap(cmap)
caxis(dlim)
set(get(hcb,'Ylabel'),'string',['Measured $',lbstr,'_{',psstr,'}$'],...
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

