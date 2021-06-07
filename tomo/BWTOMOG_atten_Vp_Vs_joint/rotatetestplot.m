typelst = {'PZmix','PZ_direct_num','P_mag_6.0_10.0','Z_mag_5.5_10.0'};
issave = true;
figure(1)
set(gcf,'position',[500 500 500 300])
for i = 1:length(typelst)
    clf;
    filename = ['rotatetest/',['rotate_test_',typelst{i},'.mat']];
    load(filename);
    plot(rotate_test.angle,rotate_test.resid,'linewidth',2);
    ylabel('$Resid$','interpreter','latex','Fontsize',12,'FontWeight','bold')
    xlabel('$Angle(deg)$','interpreter','latex','Fontsize',12,'FontWeight','bold')
    xlim([-45,0])
    grid()
    title(typelst{i},'interpreter','latex','Fontsize',12,'FontWeight','bold')
    if issave
        output = ['/media/lun/easystore/tmp/poster/rotatetest/',['rotatetest_',typelst{i}]];
        print(gcf,'-dpng',output);
    end
end