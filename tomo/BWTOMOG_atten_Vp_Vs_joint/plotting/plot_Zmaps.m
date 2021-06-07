function plot_Zmaps(plot_model,par,data,ifsave)
% plot_Zmaps(plot_model,par,data,ifsave)
HQmin = 0.3; % minimum HQ to plot

def_z = 4; % default horizontal layer to show

topo_exaggerate = 10; %3.5

%% PLOT VERTICAL SLICES

%% Velocity or Q parms
if par.t_ts==1
cmp = flipud(jet);
%cbounds = [-2.5 2.5];
cbounds = [-4 4];
elseif par.t_ts==2
cmp = parula; 
cbounds = 200*[-1 1];
end

%% Starts
grd = plot_model;
wpos=[1 500 600 530]; % window position
typestr = {'dV','dq'};

%% Bounds
min_lon = par.plot_lonlims(1); %min(mod2.lon(ind_z))+2.5;
max_lon = par.plot_lonlims(2); %max(mod2.lon(ind_z))-1.5;
min_lat = par.plot_latlims(1); %min(mod2.lat(ind_z))+2.5;
max_lat = par.plot_latlims(2); %max(mod2.lat(ind_z))-3;

%% LOAD DATA
% coast=load('~/Documents/MATLAB/CASC_atten/mapdata/m_casccoast.mat');
try
[topoX,topoY,topoZ] = grdread2('topo.grd'); % load topo grid
[topoX,topoY] = meshgrid(double(topoX),double(topoY)); topoZ = double(topoZ);
[graX,graY,graZ] = grdread2('tmp1.grd');
[graX,graY] = meshgrid(double(graX),double(graY)); graZ = double(graZ);
catch
    warning('No ETOPO data found');
end
%\jdf = dlmread('ridge_xy'); % load ridge
%fzs = dlmread('transforms_xy'); % load transforms & fracture zones
%[Vnam,Vlon,Vlat,~,~] = textread('volcanoes_complete','%s %f %f %s %s','headerlines',1);


%% GEOMETRY 
zz = unique(grd.z);
yy = [ min_lat:0.1:max_lat ]';
xx = min_lon:0.1:max_lon;


%% plot horiz slice
figure(14), clf, set(gcf,'position',[1700 200 855 1100]), hold on

% mask poor HQ
valz(grd.hq(:,:,def_z)<HQmin) = NaN;
contourf(grd.ln(:,:,def_z),grd.lt(:,:,def_z),100*grd.val(:,:,def_z),160,'edgecolor','none');
% contour HQ
[~,h] = contour(grd.ln(:,:,def_z),grd.lt(:,:,def_z),grd.hq(:,:,def_z),[HQmin:0.1:0.9]);
set(h,'LineColor','k','LineStyle','--','LineWidth',0.7);

shading flat
daspect([ 1 cosd(abs((min_lat+max_lat)/2)) 1])
colormap(cmp); caxis(cbounds); colorbar;
text(min_lon + 0.5, min_lat+0.5,[num2str(zz(def_z)),' km'],'FontSize',22,'FontWeight','bold')
% geoshow(coast.ncst(:,2), coast.ncst(:,1), 'Color', 'black','linewidth',2)
% geoshow(jdf(:,2), jdf(:,1), 'Color', 'black','linewidth',1)
% geoshow(fzs(:,2), fzs(:,1),'Linestyle','--', 'Color', 'black','linewidth',1)

% Stations
hs = plot(par.stn.lon,par.stn.lat,'v');
% set(hs,'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','c')
set(hs,'MarkerEdgeColor','k','LineWidth',.2,'MarkerSize',10,'MarkerFaceColor','none')
    
set(gca,'Layer','Top','FontSize',16)
xt=[min_lon:2:max_lon]; yt= [min_lat:2:max_lat];
set(gca,'XTick',xt,'YTick',yt) 
set(gcf,'PaperPositionMode','auto');  
xlim([min_lon max_lon])
ylim([min_lat max_lat])
box on


for ixy = 1:size(par.zsectxy,3)
    sxy = par.zsectxy(:,:,ixy);
    istr = char(64+[-1,0]+2*(ixy));
    
    figure(14); hold on
    h1 = line(sxy(:,1),sxy(:,2)); 
    set(h1,'LineWidth',2.5,'LineStyle','--','color','k')
    if diff(sxy,2) > 0 % N-S section
    text(sxy(1,1)-0.1,sxy(1,2)-0.02,istr(1),'HorizontalAlignment','center','FontWeight','bold','FontSize',22,'interpreter','latex')
    text(sxy(2,1)-0.1,sxy(2,2)+0.02,istr(2),'HorizontalAlignment','center','FontWeight','bold','FontSize',22,'interpreter','latex')
    else
    text(sxy(1,1)-0.05,sxy(1,2),istr(1),'HorizontalAlignment','right','FontWeight','bold','FontSize',22,'interpreter','latex')
    text(sxy(2,1)+0.05,sxy(2,2),istr(2),'HorizontalAlignment','left','FontWeight','bold','FontSize',22,'interpreter','latex')
    end


    % plot section
    % conditions are to get nice looking plot contextualised geographically
    figure(14+ixy), clf, set(gcf,'position',wpos + [700 0 0 0]),
    ax1 = axes('position',[0.13,0.1,0.73,0.5]);
    %axc = axes('position',[0.11,0.1,0.7,0.5]);
    ax2 = axes('position',[0.13,0.63,0.73,0.15]);
    ax3 = axes('position',[0.13,0.81,0.73,0.15]);
    ax4 = axes('position',[0.13,0.81,0.73,0.15],'YAxisLocation','right');
    
    %% get model along the section
    axes(ax4)
    %ax4.Visible = 'off';
    ax4.XTick = []; ax4.YTick = [];
    ylabel('$\leftarrow FAST\ SLOW\rightarrow$','interpreter','latex','Fontsize',8,'FontWeight','bold')
    
    [valz,hqz,slt,sln,ss,Dxy,smb,Xxy] = get_depth_section(sxy(1,:),sxy(2,:),grd);
    axes(ax3)
    
    ss = linspace(0,Xxy,length(ss));
    ss = ss-mean(ss);
    
    zthre = [200 300];
    avgv = zeros(length(ss),1);
    for k = 1:length(ss)
        tmp = valz(:,k);
        avgv(k) = sum(tmp(zz<=zthre(2)&zz>=zthre(1)))/sum(zz<=zthre(2)&zz>=zthre(1));
    end
    try
    plot(ss,avgv*100,'k','LineWidth',2);
    ylim([-max(abs(avgv))*100,max(abs(avgv))*100])
    catch
    plot(ss,zeros(length(ss),1),'k','LineWidth',2);
    end
    xlim([min(ss) max(ss)])
    set(ax3,'Xtick',[-200:50:200],'Xticklabel',[],'Layer','Top','FontSize',8,'ydir','reverse')
    ylabel('$avg\ \delta V\%$','interpreter','latex','Fontsize',12,'FontWeight','bold')
    
   
    set(ax3,'Xtick',[-200:50:200],'Xticklabel',[],'Layer','Top','FontSize',8,'ydir','reverse')
    ylabel('$avg\ \delta V\%$','interpreter','latex','Fontsize',12,'FontWeight','bold')
    % % mask out low hq
    valz(hqz<HQmin) = NaN;
    % hqz(zz<40,:)=NaN;
    
    axes(ax1); hold on;
    contourf(ss,zz,100*valz,160,'edgecolor','none');
    set(gca,'YDir','reverse','FontSize',12);
    %axis equal

    shading flat
    set(gca,'Layer','Top','FontSize',8,...
        'YTick',[zz(1),50:50:zz(end)]','YTickLabel',num2str([zz(1),50:50:zz(end)]'))
    ylabel('$Depth(km)$','FontSize',12,'FontWeight','bold','interpreter','latex')
    xlabel('$Distance\ along\ profile(km)$','FontSize',12,'FontWeight','bold','interpreter','latex')
    colormap(cmp); caxis(cbounds);
    h1 = colorbar(ax1,'Position',[.88,.1,.03,.5],'AxisLocation','in');
    set(h1,'FontSize',8,'YAxisLocation','right','Ytick',linspace(cbounds(1),cbounds(2),5))
    str = ['$\delta V \,\, \%$'];
    ylabel(h1,str,'FontSize',12, 'Interpreter','Latex')
    % contour HQ
    %[~,h] = contour(ss,zz,hqz,[0.6:0.1:0.9]);
    %set(h,'LineColor','k','LineStyle','--','LineWidth',0.7);

    % % contour semb
    % [cs,h] = contour(ss,zz,smb,0.7);
    % set(h,'LineColor','b','LineStyle','--','LineWidth',2);
    

    % CHANGE SCALE AT 5 KM DEPTH
    %shifter = 5*(topo_exaggerate - 1); % from simple maths

    %% plot on topography

for i=1:1
    axes(ax2);
    [dis0 azi0] = distance(slt(1),sln(1),slt(end),sln(end));
    [slt(1),sln(1)] = latlon_from(slt(1),sln(1),azi0+180,270);
    span = 100.;
    gap = 4.;
    [lat1b,lon1b] = latlon_from(slt(1),sln(1),azi0-90,gap:gap:span);
    [lat2b,lon2b] = latlon_from(slt(1),sln(1),azi0+90,gap:gap:span);

    lat = [lat1b slt(1) lat2b];
    lon = [lon1b sln(1) lon2b];
    %lat = [slt(1)];
    %lon = [sln(1)];
    Profile_width = 1;
    Profile_range = (0:Profile_width:deg2km(dis0)+540);
    
    avg_topo = zeros(1,length(Profile_range));
    avg_gra = zeros(1,length(Profile_range));
    for ip = 1:size(lat,2)
	[Profile_lat, Profile_lon] = latlon_from(lat(ip),lon(ip),azi0,Profile_range);
	avg_topo = avg_topo + interp2(topoX,topoY,topoZ,Profile_lon,Profile_lat)/1000;
	avg_gra = avg_gra + interp2(graX,graY,graZ,Profile_lon,Profile_lat);
    end
    avg_topo = avg_topo/size(lat,2);
    avg_gra = avg_gra/size(lat,2);
    %topo = interp2(topoX,topoY,topoZ,sln,slt)/1000;
    %gra = interp2(graX,graY,graZ,sln,slt);
    %topo = -topo./(1000./topo_exaggerate);
        % CHANGE SCALE AT 5 KM DEPTH
    %topo = topo-shifter;
    dis_axis = linspace(-mean(Profile_range),mean(Profile_range),length(Profile_range));
    [Profile_lat, Profile_lon] = latlon_from(slt(1),sln(1),azi0,Profile_range);
    
    save('profile.mat','avg_topo','avg_gra','dis_axis','Profile_lat','Profile_lon');
    
    [AX,H1,H2] = plotyy(dis_axis,avg_topo,dis_axis,avg_gra);
    set(H1,'LineWidth',2);
    set(H2,'LineWidth',2);
    set(get(AX(1),'Ylabel'),'String','$Topo(km)$');
    set(get(AX(2),'Ylabel'),'String','$Gravity(mgal)$');
    set(AX,'Xlim',[min(ss) max(ss)],'Xtick',[-200:50:200],'Xticklabel',[],'FontSize',8)
    set(AX(1),'Ylim',[-0.1 0.1],'Ytick',-0.1:0.1:0.1)
    set(AX(2),'Ylim',[-5,4.6],'Ytick',[-4:2:4])
    %y1 = floor(min(abs(min(avg_topo)),abs(max(avg_topo)))*10)/10.;
    %set(AX(1),'Ylim',[min(avg_topo) max(avg_topo)],'Ytick',[-y1 0 y1])
    %y2 = floor(min(abs(min(avg_gra)),abs(max(avg_gra)))*10)/10.;
    %set(AX(2),'Ylim',[min(avg_gra) max(avg_gra)],'Ytick',[-y2 0 y2])
    ylabel(AX(1),'$Topo(km)$','interpreter','latex','Fontsize',12,'FontWeight','bold')
    ylabel(AX(2),'$Gra(mgal)$','interpreter','latex','Fontsize',12,'FontWeight','bold')
    %ht = plot(ss,topo/1000,'k','LineWidth',3);
    %set(ax2,'Xtick',[-200:100:200],'Xticklabel',[],'Layer','Top','FontSize',8)
    %ylabel('$Topo(km)$','interpreter','latex','Fontsize',12,'FontWeight','bold')

    % plot water
    %water.x = [ss,fliplr(ss)]';
    % CHANGE SCALE AT 5 KM DEPTH (otherwise would be water.y = [zeros(size(topo)); flipud(topo)]; water.y(water.y<0)=0;
    %water.y = [-shifter*ones(size(topo)); flipud(topo)]; water.y(water.y<-shifter)=-shifter;
    %hw = patch(water.x,water.y,'b');
    %set(hw,'EdgeAlpha',0,'FaceAlpha',0.5)
end
%catch
%    warning('No ETOPO data found');
%end
%{
    %% plot on slab
    %[slbs,slbz] = slab_on_section(sxy);
    %hsl = plot(slbs-(Xxy/2),-slbz,'--k','LineWidth',3);

    
    %% plot volcanoes within 50 km
%     % calc. distance of stations to line, which are in bounds
%     Dev = dist2line(sxy(1,:),sxy(2,:),[Vlon,Vlat])*Xxy/norm(diff(sxy));
%     indx = find(abs(Dev)<50);
%     % calc. projection of each stations along line
%     Sev = ([Vlon,Vlat]-ones(length(Vlon),1)*sxy(1,:))*diff(sxy)'*Xxy/(norm(diff(sxy))^2);
%     % only include stations within end-bounds of line
%     indx1 = find(Sev<=Xxy & Sev>=0);
%     indx = intersect(indx,indx1);
%     hs = plot(Sev(indx)-(Xxy/2),interp1(ss,topo,Sev(indx)-(Xxy/2)),'^');
%     set(hs,'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',12,'MarkerFaceColor','r')
%}
    %% PLOT STATIONS within 50 km
    axes(ax1)
    %{
    % calc. distance of stations to line, which are in bounds
    Dev = dist2line(sxy(1,:),sxy(2,:),[data.stn.lon,data.stn.lat])*Xxy/norm(diff(sxy));
    indx = find(abs(Dev)<50);
    % calc. projection of each stations along line
    Sev = ([data.stn.lon,data.stn.lat]-ones(length(data.stn.lon),1)*sxy(1,:))*diff(sxy)'*Xxy/(norm(diff(sxy))^2);
    % only include stations within end-bounds of line
    indx1 = find(Sev<=Xxy & Sev>=0);
    indx = intersect(indx,indx1);
    %hs = plot(Sev(indx)-(Xxy/2),-data.stn.elv(indx)*topo_exaggerate-shifter,'v');
    hs = plot(Sev(indx)-(Xxy/2),4.5,'v');
    set(hs,'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','c')
    uistack(hs,'top');
    try
    uistack(hw,'bottom');
    uistack(ht,'bottom')
    end
    %}
    
    
    %% plot box and axes
    %set(gca,'Layer','Top','FontSize',16,...
    %    'YTick',[-shifter,5,50:50:par.max_z]','YTickLabel',num2str([0,5,50:50:par.max_z]'))
    %ylim([-shifter - 5*topo_exaggerate,zz(end)])
    ylim([zz(1) zz(end)])
    xlim([min(ss), max(ss)])
    box on

    % plot section ends
    text(min(ss)+10,zz(1)+10,istr(1),'HorizontalAlignment','center','FontWeight','bold','FontSize',16)
    text(max(ss)-10,zz(1)+10,istr(2),'HorizontalAlignment','center','FontWeight','bold','FontSize',16)

    axes(ax1)
    
    %% SAVE
    set(gcf, 'Color', [1 1 1])
    if ifsave
        save2pdf(14+ixy,strcat(typestr{par.t_ts},'_Zmap',num2str(ixy)),par.figdir);
    end

end % loop on sections

%% Save hmap
    if ifsave
        save2pdf(14,strcat(typestr{par.t_ts},'_Zmap_hplot'),par.figdir);
    end

end

function [slbs,slbz] = slab_on_section(sxy)
    Xxy = distance_km(sxy(1,2),sxy(1,1),sxy(2,2),sxy(2,1));
    [ slab_contrs ] = McCrory_slab_depth;
    Nc = length(slab_contrs);
    slbs = nan(Nc,1);
    slbz = nan(Nc,1);
    for ic = 1:Nc
        dst = dist2line(sxy(1,:),sxy(2,:),[slab_contrs(ic).lon,slab_contrs(ic).lat]);
        if ~any(dst<0), continue; end
        if ~any(dst>0), continue; end
        ind = [find(dst==min(dst(dst>0))),find(dst==max(dst(dst<0)))];        
        ss = ([slab_contrs(ic).lon(ind),slab_contrs(ic).lat(ind)]-ones(2,1)*sxy(1,:))*diff(sxy)'*Xxy/(norm(diff(sxy))^2);
        slbs(ic) = mean(ss);
        slbz(ic) = slab_contrs(ic).depth;
    end
    
    % fill in any nans w/ linear interp
    nnan = ~isnan(slbz);
    slbs = interp1(slbz(nnan),slbs(nnan),[slab_contrs.depth]);
    slbz = interp1(slbz(nnan),slbz(nnan),[slab_contrs.depth]);
end

