function plot_Hmaps(plot_model,par,opt,saveopt)
% opt is an option describing the input model
%  opt == 1  means real result
%  opt == 2  means synth in
%  opt == 3  means synth out
%  opt == 4  means plot bootstrap uncertainties
%
% saveopt is an option to save (1) or not (0)
if nargin < 3
    opt = 1;
end
if nargin < 4
    saveopt = 0;
end

if par.t_ts == 1
    cmp = flipud(jet);
    valstr = 'V';
    %clims = [-2.5 2.5];
    clims = [-4 4];
elseif par.t_ts == 2
    cmp = parula;
    valstr = 'q';
    clims = 100*[-1 1];
end

if par.PS == 3 && par.Rdvpdvs==0
    mdim = 2; phstr = {'P','S'}; vmul = [1,2];
elseif par.PS == 3 && par.Rdvpdvs~=0
    mdim = 1; phstr = {'S'}; vmul = [2];
elseif par.PS == 1
    mdim = 1; phstr = {'P'}; vmul = [1];
elseif par.PS == 2
    mdim = 1; phstr = {'S'}; vmul = [2];
end

% min_lon = -132; %min(mod2.lon(ind_z))+2.5;
% max_lon = -119; %max(mod2.lon(ind_z))-1.5;
% min_lat = 39; %min(mod2.lat(ind_z))+2.5;
% max_lat = 51; %max(mod2.lat(ind_z))-3;

HQmin = 0.3;
%zz = [100,120,140,160,180,200,220,240,260];
zz = [6,40,80,110,140,180,220,260,300];
%zz = 40:20:200;
nz = length(zz);

try
[graX,graY,gra] = grdread2('plotting/tmp1.grd'); % load topo grid
[graX,graY] = meshgrid(double(graX),double(graY)); gra = double(gra);
catch
    warning('No ETOPO data found');
end


%% load some data
% mapdata = '/Users/zeilon/Documents/MATLAB/CASC_atten/mapdata/';
% coast = load([mapdata,'m_casccoast.mat']); % load coastline
% jdf = dlmread([mapdata,'ridge_xy']); % load ridge
% fzs = dlmread([mapdata,'transforms_xy']); % load transforms & fracture zones


%% ===========================  PLOT VELOCITY  ===========================
%% ===========================  PLOT VELOCITY  ===========================

% loop over model dimensions (P and S) if applicable
for ip = 1:mdim
    
clims = vmul(ip)*clims;

fign = 57 + 10*opt + 100*(ip-1);

figure(fign); clf

set(gcf,'Position',[0,280-20*opt,900,900])
%[nxsub,nysub] = nsubplots(nz);
%pos = fun_mm_subplot_pos(3,3,0.95,0.88);
for iz = 1:nz
    %subplot(nysub,nxsub,iz)
    ax{iz} = axes('position',[0.07+0.26*mod(iz-1,3),0.68-0.22*floor((iz-1)/3),0.24,0.22]);
    box;
	hold on
    
    yy = par.plot_latlims(1):0.1:par.plot_latlims(2);
    xx = par.plot_lonlims(1):0.1:par.plot_lonlims(2);
    z = zz(iz);
    [xmesh,ymesh,zmesh] = meshgrid(xx,yy,z);

    if opt==2
	for i=1:par.nz-1
	    if (par.zz(i)<=zz(iz))&&(par.zz(i+1)>zz(iz))
		j = i;
	        break
	    end
	end
	VAL = griddata(plot_model.ln(:,:,j),plot_model.lt(:,:,j),100*plot_model.val(:,:,j,ip),xx,yy');
        cbounds = clims;

    elseif opt~=4

        VAL = griddata(plot_model.ln,plot_model.lt,plot_model.z,100*plot_model.val(:,:,:,ip),xmesh,ymesh,zmesh);
        cbounds = clims;
        
    elseif opt==4
       
        VAL = griddata(plot_model.ln,plot_model.lt,plot_model.z,100*plot_model.sv,xmesh,ymesh,zmesh);
        cbounds = [0 0.1];
        cmp = colormap(cmap_makecustom([0.8 0.8 0.1],[0.1 0.8 0.8],0));
        
    end
%     smb = griddata(plot_model.ln(:,:,iz),plot_model.lt(:,:,iz),plot_model.semb_v(:,:,iz),xx,yy);
    %hq = griddata(plot_model.ln(:,:,iz),plot_model.lt(:,:,iz),plot_model.hq(:,:,iz,ip),xx,yy);
    hq = griddata(plot_model.ln,plot_model.lt,plot_model.z,plot_model.hq(:,:,:,ip),xmesh,ymesh,zmesh);
        
    VAL(hq<HQmin) = nan;
    
    contourf(xx,yy,VAL,80,'edgecolor','none');
    [~,gracontour] = contour(graX,graY,gra,[-20:6:15]);
    %set(gracontour,'LineColor','k','LineStyle','-','LineWidth',0.7);
%     contour(xx,yy,smb,[0.7:0.1:1],'--r','Linewidth',1.5)
%     contour(xx,yy,smb.*hq,[0.5:0.1:1],'--b','Linewidth',1.5)
    hs = plot(par.stn.lon,par.stn.lat,'v');
    set(hs,'MarkerEdgeColor','k','LineWidth',.2,'MarkerSize',5,'MarkerFaceColor','none')
    shading flat
    axis( [par.plot_lonlims par.plot_latlims] );
    daspect([ 1 cosd(abs(mean(par.plot_latlims))) 1])
    colormap(cmp)
    caxis(cbounds)
    
%     geoshow(coast.ncst(:,2), coast.ncst(:,1), 'Color', 'black','linewidth',2)
%     geoshow(jdf(:,2), jdf(:,1), 'Color', 'black','linewidth',1)
%     geoshow(fzs(:,2), fzs(:,1),'Linestyle','--', 'Color', 'black','linewidth',1)
    
%     xlabel('Longitude'); 
%     ylabel('Latitude');
    xt = par.plot_lonlims(2) - 0.25*(par.plot_lonlims(2) - par.plot_lonlims(1));
    yt = par.plot_latlims(2) - 0.1*(par.plot_latlims(2) - par.plot_latlims(1));
    text(xt,yt,sprintf('%.0fkm',zz(iz)),'FontSize',10,'FontWeight','bold');
    if mod(iz-1,3) ~= 0
        set(ax{iz},'yTick',[])
    else
        ylabel('$Latitude[^{\circ}]$','FontWeight','bold','FontSize',12,'interpreter','latex')
    end
    if iz<= 6
        set(ax{iz},'xTick',[])
    else
        xlabel('$Longitude[^{\circ}]$','FontWeight','bold','FontSize',12,'interpreter','latex')
	set(ax{iz},'xTick',-136:-131)	
    end
    %title(sprintf('Depth slice %.0f km',zz(iz)),'FontSize',14)
end

%% scale


%subplot(nxsub,nysub,iz+1)
%axes('position',[0.93,0.1,0.05,0.8]);
%set(gca,'Visible','off')
if opt~=4
    str = ['$\delta ',valstr,' \,\, \%$']; 
    colormap(cmp); 
    caxis(clims)
elseif opt==4
    str = '%\sigma_{boot}$'; 
    colormap(cmap_makecustom([0.8 0.8 0.1],[0.1 0.8 0.8])); 
    caxis([0 0.1]);
end
hc = colorbar('EastOutside');
set(get(hc,'YLabel'),'String',str,'FontSize',20,'FontWeight','bold','interpreter','latex')
set(hc,'FontSize',12,'FontWeight','bold','Position',[0.84 0.25 0.02 0.64])

%% title
%title_custom( sprintf('V%s anomaly',phstr{ip}),0.97,0.5)

%% save
if saveopt
    suff = {'all','all_synin','all_synout','booterrs'};
    pref = ['d',valstr];
    if par.wtdata == 0, wtstr = '_nowt'; else wtstr = '';  end
   
    fprintf('>  Saving figure %s_%s%s... \n',pref,suff{opt},wtstr);
    
    ostr = sprintf('%s/%s_%s%s.eps',par.figdir,pref,suff{opt},wtstr);
    print(ostr,'-depsc','-r600');

end

pause(0.01)

end

%% subfunctions
