fold = 'sqztest';
tmp = dir(fold);
issave = true;
figure(1)
set(gcf,'position',[500 500 500 300])
for i = 3:length(tmp)
    clf;
    titletxt = replace(tmp(i).name,'.mat','');
    filename = fullfile(fold,tmp(i).name);
    load(filename);
    [AX,H1,H2] = plotyy(squeezetest.zsqz,squeezetest.wvr./squeezetest.nmodel',squeezetest.zsqz,squeezetest.m2norm);
    set(AX,'position',[0.15,0.1,0.72,0.8])
    set(H1,'LineWidth',2);
    set(H2,'LineWidth',2);
    set(AX,'Xlim',[min(squeezetest.zsqz) max(squeezetest.zsqz)],'FontSize',8)
    set(AX(1),'ytick',ceil(min(squeezetest.wvr./squeezetest.nmodel')*500)/500:floor((max(squeezetest.wvr./squeezetest.nmodel')-min(squeezetest.wvr./squeezetest.nmodel'))*500)/2500:max(squeezetest.wvr./squeezetest.nmodel'))
    %set(AX(1),'ytick',ceil(min(squeezetest.wvr)):floor((max(squeezetest.wvr)-min(squeezetest.wvr))/5):max(squeezetest.wvr))
    ylabel(AX(1),'$Variance\ Reduction(\%)$','interpreter','latex','Fontsize',12,'FontWeight','bold')
    ylim(AX(1),[min(squeezetest.wvr./squeezetest.nmodel') max(squeezetest.wvr./squeezetest.nmodel')])
    set(AX(2),'ytick',ceil(min(squeezetest.m2norm)*500)/500:floor((max(squeezetest.m2norm)-min(squeezetest.m2norm))*500)/2500:max(squeezetest.m2norm))
    ylabel(AX(2),'$m2\ norm$','interpreter','latex','Fontsize',12,'FontWeight','bold')
    ylim(AX(2),[min(squeezetest.m2norm) max(squeezetest.m2norm)])
    xlabel('$zsqz(km)$','interpreter','latex','Fontsize',12,'FontWeight','bold')
    %grid()
    title(titletxt,'interpreter','latex','Fontsize',12,'FontWeight','bold')
    if issave
        output = ['/media/lun/easystore/tmp/poster/sqztest/sqzcurve_normalizedvr/',titletxt,'.png'];
        print(gcf,'-dpng',output);
    end
end
    

%{
typelst = {'PZmix','Z_mag_5.5_10.0','PZ_direct_num','P_mag_6.0_10.0'};
is2D = false;
issave = true;
figure(1)
set(gcf,'position',[500 500 500 300])
for i = 1:length(typelst)
    clf;
    if ~is2D
        filename = ['sqztest/',['sqztest_',typelst{i},'.mat']];
        titletxt = typelst{i};
    else
        filename = ['sqztest/',['sqztest_',typelst{i},'_2D.mat']];
        titletxt = [typelst{i},'_2D'];
    end
    load(filename);
    [AX,H1,H2] = plotyy(squeezetest.zsqz,squeezetest.wvr,squeezetest.zsqz,squeezetest.m2norm);
    set(H1,'LineWidth',2);
    set(H2,'LineWidth',2);
    set(AX,'Xlim',[min(squeezetest.zsqz) max(squeezetest.zsqz)],'FontSize',8)
    ylabel(AX(1),'$Variance\ Reduction(\%)$','interpreter','latex','Fontsize',12,'FontWeight','bold')
    ylabel(AX(2),'$m2\ norm$','interpreter','latex','Fontsize',12,'FontWeight','bold')
    xlabel('$zsqz(km)$','interpreter','latex','Fontsize',12,'FontWeight','bold')
    grid()
    title(titletxt,'interpreter','latex','Fontsize',12,'FontWeight','bold')
    if issave
        output = ['/media/lun/easystore/tmp/IM2.22/sqztest/',titletxt];
        print(gcf,'-dpng',output);
    end
end
%}