function figh = plot_sta_static(par,data,model,inmodel)
% figh = plot_sta_static(par,data,model,inmodel)

% defaults
isin = 0;
nph = 1;

if nargin > 3 && ~isempty(inmodel)
    isin = 1;
end 
if par.PS == 3 
    nph = 2; psstr = {'P','S'};
elseif par.PS == 2 
    psstr = {'S'};
elseif par.PS == 1
    psstr = {'P'};
end

figh = figure(88); clf;
set(figh,'pos',[30 20 580*nph 390*(1+isin)])

for iph = 1:nph
    x0 = 0.1/nph + (iph-1)*0.5;
    dx = 0.85/nph;
    y0 = 0.12/(1+isin);
    dy = 0.82/(1+isin);
    ax(iph,1) = axes('pos',[x0,y0,dx,dy]);

    scatter(ax(iph,1),data.stn.lon,data.stn.lat,50,...
        -model.sstatic([1:data.stn.nstas] + data.stn.nstas*(iph-1)),'filled')
    title(ax(iph,1),sprintf('Final %s S-static',psstr{iph}))
    
    if isin
        ax(iph,2) = axes('pos',[x0,y0 + 1.15*dy,dx,dy]);
        scatter(ax(iph,2),data.stn.lon,data.stn.lat,50,...
            inmodel.sstatic([1:data.stn.nstas] + data.stn.nstas*(iph-1)),'filled'), 
        title(ax(iph,2),sprintf('Input %s S-static',psstr{iph}))
    end
    xlabel('lon')
    ylabel('lat')
end

clim = [-0.5,0.5];
cmp = jet;
for ix = 1:numel(ax)
    axis(ax(ix),[par.plot_lonlims par.plot_latlims])
    colormap(cmp)
    caxis(clim)
    h = colorbar;
    h.Label.String = 'dT[s]'
    %set(h, 'ylim', [-0.5 0.5])
end

end

