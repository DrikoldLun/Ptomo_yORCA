%% %% %% 
% clear all
addpath('matguts','plotting','function','seizmo');

saveopt = 1;
smoothtype = 'Ltest_flatness';
%smoothtype = 'Ltest_smoothness';
savepath = ['/media/lun/easystore/tmp/poster/Ltest/',smoothtype];
%ofile = 'Ltest_dT_PZmix.mat'; %0.6 0.5 0.6 0.6
%ofile = 'Ltest_dT_PZ_direct_num.mat'; %0.6 0.6 0.55 0.55
%ofile = 'Ltest_dT_Z_mag_5.5_10.0.mat'; %0.4 0.38 0.4 0.4
ofile = 'Ltest_dT_P_mag_6.0_10.0.mat'; %2 1.8 2.1 1.7
type = replace(replace(ofile,'.mat',''),'Ltest_dT_','');
ofile = fullfile(smoothtype,ofile);
ofig = 'Lcurv_dT_PZ';
figN = 6;
resid2rough = 2.1;
load(ofile)
Ltest.vr = Ltest.wvr;


Nd = length(Ltest.damp);
Ns = length(Ltest.smooth);
cls_d = colour_get(1:Nd,Nd,1,flipud(autumn));
cls_s = colour_get(1:Ns,Ns,1,flipud(winter)); %parula
close all

%% Calc norms
%norm = Ltest.norm + Ltest.norm_estt + Ltest.norm_sstt;
%norm = Ltest.smth;
%norm = (Ltest.smth./max(Ltest.smth) + Ltest.norm./max(Ltest.norm))*100;
norm = (Ltest.smth./max(Ltest.smth(:)) + 0.1*Ltest.norm./max(Ltest.norm(:)))*100;
%norm = Ltest.smth + 0.2*Ltest.norm;

%% plot
figure(33), clf
set(gcf,'Position',[30 330 700 500])
ax1 = axes;
ax2 = axes;
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
set([ax1,ax2],'Position',[.1 .1 .75 .75]);
%ax1 = axes('position',[0.13 0.19 0.8 0.75]);
axes(ax1), hold on

for iss = 1:Ns
% scatter(Ltest.vr(:,iss),norm_all(:,iss),50,cls_d,'filled')
% scatter(Ltest.vr_tt(:,iss),norm_tt(:,iss),50,cls_d,'filled')
% scatter(Ltest.vr_dt(:,iss),norm_dT(:,iss),50,cls_d,'filled')  
for idd = 1:Nd
    plot(Ltest.vr(idd,iss),norm(idd,iss),'o','MarkerEdgeColor','none','MarkerSize',10,'MarkerFaceColor',cls_d(idd,:))
%     plot(Ltest.vr_tt(idd,iss),norm_tt(idd,iss),'o','MarkerEdgeColor','none','MarkerSize',5,'MarkerFaceColor',cls_d(idd,:))
%     plot(Ltest.vr_dt(idd,iss),norm_dT(idd,iss),'o','MarkerEdgeColor','none','MarkerSize',5,'MarkerFaceColor',cls_d(idd,:))
end

end
for iss = 1:Ns
    plot(Ltest.vr(:,iss),norm(:,iss),'-','MarkerSize',4,'Linewidth',1.,'color',cls_s(iss,:));
end

%plot(Ltest.resid(Ltest.damp==3,Ltest.smooth==3),norm(Ltest.damp==3,Ltest.smooth==3),...
%    'ok','MarkerSize',15,'MarkerFaceColor','k')

%% axis etc.
%set(ax1,'Fontsize',16,'xdir','reverse','LineWidth',1.5,'XTick',[0:10:100],'YTick',[0:2:20])
set(ax1,'Fontsize',16,'xdir','reverse','LineWidth',1.5)
xlabel('Variance reduction(\%)',  'FontSize',12,'Interpreter','Latex')
%ylabel('Model norm (s)',          'FontSize',23,'Interpreter','Latex')
%ylabel('Normalized roughness+model_length(s)',          'FontSize',12,'Interpreter','Latex')
ylabel('X',          'FontSize',12,'Interpreter','Latex')
%set(get(ax1,'Ylabel'),'Position',[87.8 3.5 1])
%axis([10 100 0 11])

% text(101,-18,['\textbf{Figure ',num2str(figN),'}'],'Fontsize',25,'interpreter','latex')


%% legend + annotations
vrmax = max(Ltest.vr(:));vrmin = min(Ltest.vr(:));
%vrmax = 100; vrmin = 10;
vrspan = vrmax-vrmin;
normmax = max(norm(:));normmin = min(norm(:));
%normmax = 11; normmin = 0;
normspan = normmax-normmin;
axis([vrmin vrmax normmin normmax])

text(0.95*vrspan+vrmin,0.97*normspan+normmin,'$\leftarrow$ less damped','interpreter','latex','fontsize',15,'rotation',-45)
text(0.26*vrspan+vrmin,0.05*normspan+normmin,'more damped $\rightarrow$ ','interpreter','latex','fontsize',15,'rotation',-20)
%text(97,8,'$\gamma = 3$','interpreter','latex','fontsize',20)
%text(96.7,14,'$\epsilon=3$','interpreter','latex','fontsize',20)
%arrow([92,17],[87.2,34],'length',10)

%% key
keyloc = [0.15*vrspan+vrmin,0.75*normspan+normmin];
keysiz = [0.2*vrspan,0.35*normspan];

kle = keyloc(1) - 0.5*keysiz(1); kri = keyloc(1) + 0.5*keysiz(1);
kbo = keyloc(2) - 0.5*keysiz(2); kto = keyloc(2) + 0.5*keysiz(2);
kw = keysiz(1); kh = keysiz(2);
%{
hold on
patch([kle kle kri kri kle],[kbo kto kto kbo kbo],...
      'w','LineWidth',2)
text(keyloc(1),kto-0.03*kh,'smoothing ($\gamma$)',...
    'fontsize',12,'fontweight','bold','verticalalignment','top','horizontalalignment','center','interpreter','latex');

for iss = 1:4:21
plot(kle + [0.3*kw,0.85*kw],kbo + (iss/(Ns+3))*kh*0.85 + 0.03*kh + [0 0],'LineWidth',2,'color',cls_s(iss,:))
text(kle + 0.2*kw,kbo + (iss/(Ns+3))*kh*0.85 + 0.03*kh ,num2str(Ltest.smooth(iss)),'verticalalignment','middle','FontSize',6,'fontweight','bold')
end
%}


% [.1 .1 .75 .75]


 
% return
% scatter(kle + 0.3*kw,kbo + 0.61*kh,sqrt(1),  'k','filled','MarkerEdgeColor','k');
% scatter(kle + 0.3*kw,kbo + 0.435*kh,sqrt(10), 'k','filled','MarkerEdgeColor','k');
% scatter(kle + 0.3*kw,kbo + 0.18*kh,sqrt(50),'k','filled','MarkerEdgeColor','k');
% 
% text(kle + 0.62*kw,kbo + 0.61*kh,'1', 'fontsize',12,'fontweight','bold','verticalalignment','middle');
% text(kle + 0.62*kw,kbo + 0.435*kh,'10','fontsize',12,'fontweight','bold','verticalalignment','middle');
% text(kle + 0.62*kw,kbo + 0.18*kh,'50','fontsize',12,'fontweight','bold','verticalalignment','middle');




%% penalty function
penalty = norm + resid2rough*(100-Ltest.vr);
%penalty = (Ltest.norm) + 0.5*(100-Ltest.vr);

% subplot(2,1,2)
figure(34), clf
set(gcf,'Position',[100 600 600 450])
%h = contourf(Ltest.smooth,Ltest.damp,penalty,100);
imagesc(Ltest.smooth,Ltest.damp,penalty);
set(gca,'ydir','normal')
shading flat
%set(gca,'XTick',1:length(Ltest.smooth),'XTickLabel',Ltest.smooth,...
%        'YTick',1:4:length(Ltest.damp),'YTickLabel',Ltest.damp(1:4:end))
%set(gca,'Fontsize',18,'yscale','log')
xlabel('Smoothing ($\gamma$)',  'FontSize',12,'Interpreter','Latex')
ylabel('Damping ($\epsilon$)',  'FontSize',12,'Interpreter','Latex')


[minA,x,y] = mingrid(penalty);
title(sprintf('Contour plot of penalty: (F = %.3f)\nMinimum is damp = %.1f, smooth = %.1f',...
    resid2rough,Ltest.damp(y),Ltest.smooth(x)),'FontSize',12, 'Interpreter','Latex')
%set(get(gca,'Ylabel'),'Position',[0.27 7 1])

clim = [minA, 1.15*minA];
cmp = jet;
colormap(cmp)
caxis(clim)
h = colorbar;
ylabel(h,'Penalty','FontSize',12, 'Interpreter','Latex')
%caxis([minA, 1.5*minA]);

%% plot the minimum:
% subplot(2,1,1), 
figure(33), hold on
plot(Ltest.vr(y,x),norm(y,x),'o','MarkerSize',10,'MarkerEdge','w','MarkerFace','k')
xline = linspace(0,100,10000);
yline = linspace(minA-resid2rough*100, minA+resid2rough*100-resid2rough*100, 10000);
line = plot(xline, yline, 'k-', 'linewidth', 3);
legend(line,'Penalty function','FontSize',12,'Interpreter','Latex')
set(gca,'yscale','log')
ax1.FontSize = 10;

clim = [1,Ns];
cmp = flipud(winter);
colormap(ax1,cmp)
caxis(clim)
h1 = colorbar(ax1,'Position',[.87,.1,.03,.75],'AxisLocation','in','Location','east');
set(h1, 'YAxisLocation','right')
ylabel(h1,'Smooth','FontSize',12, 'Interpreter','Latex')
h1.FontSize = 10;
ytick = [1:4:21];
for i = 1:length(ytick)
    yticklabel{i} = num2str(Ltest.smooth(ytick(i)));
end
set(h1,'YTick',[1:4:21],'Yticklabel',yticklabel)

axes(ax2), hold on
clim = [1,Nd];
cmp = flipud(autumn);
colormap(ax2,cmp)
caxis(clim)
h2 = colorbar(ax2,'Position',[.1,.87,.75,.03],'AxisLocation','in','Location','north');
h2.FontSize = 10;
set(h2, 'YAxisLocation','top')
ylabel(h2,'Damp','FontSize',12, 'Interpreter','Latex')
ytick = [1:4:21];
for i = 1:length(ytick)
    yticklabel{i} = num2str(Ltest.damp(ytick(i)));
end
set(h2,'YTick',[1:4:21],'Yticklabel',yticklabel)

if saveopt
    print(gcf,'-dpng',fullfile(savepath,sprintf('Lcurve_%s_%s.png',type,smoothtype)));
end
% subplot(2,1,2), hold on
figure(34), hold on
plot(Ltest.smooth(x),Ltest.damp(y),'o','MarkerSize',14,'MarkerEdge','w','MarkerFace','k')
if saveopt
    print(gcf,'-dpng',fullfile(savepath,sprintf('penalty_%s_%s.png',type,smoothtype)));
end

return
%{
if saveopt
    save2pdf(33,ofig,'~/Documents/MATLAB/CASC_atten/TOMOGRAPHY/figs/');
end

if saveopt
    save(ofile,'Ltest');
    save2pdf(34,ofig,'~/Documents/MATLAB/CASC_atten/TOMOGRAPHY/figs/');
end
%}


