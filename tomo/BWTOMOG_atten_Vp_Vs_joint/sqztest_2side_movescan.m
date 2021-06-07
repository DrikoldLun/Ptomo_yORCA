load('sqztest_2side/sqz2side_PZ_direct_num.mat')
zs = squeezetest(:).zs;
zd = squeezetest(:).zd;
wvr = squeezetest(:).wvr;
m2norm = squeezetest(:).m2norm;
dz = zd - zs;
%ind = find(dz==0); %1-layer
%ind = find(dz>10&dz<50); %2-layer
%ind = find(dz>40&dz<90); %3-layer
ind = find(dz>90&dz<130); %4-layer
zs = zs(ind);
zd = zd(ind);
wvr = wvr(ind);
m2norm = m2norm(ind);

figure(1)
clf;
set(gcf,'position',[500 500 500 300])
[AX,H1,H2] = plotyy(zs,wvr,zs,m2norm);
%set([H1;H2],'Marker','.','MarkerSize',24)
set(AX,'position',[0.15,0.12,0.72,0.8])
set(H1,'LineWidth',2);
set(H2,'LineWidth',2);
set(AX,'Xlim',[min(squeezetest(:).zs) max(squeezetest(:).zs)],'FontSize',8,'Xtick',unique(squeezetest(:).zs))
xlabel(AX(1),'$Z_s(km)$','interpreter','latex','Fontsize',12,'FontWeight','bold')

%set(AX(1),'ytick',ceil(min(squeezetest.wvr./squeezetest.nmodel')*500)/500:floor((max(squeezetest.wvr./squeezetest.nmodel')-min(squeezetest.wvr./squeezetest.nmodel'))*500)/2500:max(squeezetest.wvr./squeezetest.nmodel'))
%set(AX(1),'ytick',ceil(min(squeezetest.wvr)):floor((max(squeezetest.wvr)-min(squeezetest.wvr))/5):max(squeezetest.wvr))
axes(AX(1))
hold on
plot(zs,wvr,'Marker','.','MarkerSize',30,'color','blue')
set(AX(1),'ytick',linspace(65,86,5))
ylabel(AX(1),'$Variance\ Reduction(\%)$','interpreter','latex','Fontsize',12,'FontWeight','bold')
ylim(AX(1),[65 86])

axes(AX(2))
hold on
plot(zs,m2norm,'Marker','.','MarkerSize',30,'color',[0.8500, 0.3250, 0.0980])
%set(AX(2),'ytick',ceil(min(squeezetest.m2norm)*500)/500:floor((max(squeezetest.m2norm)-min(squeezetest.m2norm))*500)/2500:max(squeezetest.m2norm))
ylabel(AX(2),'$m2\ norm$','interpreter','latex','Fontsize',12,'FontWeight','bold')
ylim(AX(2),[0.225 0.345])
set(AX(2),'ytick',linspace(0.225,0.345,5))
grid()