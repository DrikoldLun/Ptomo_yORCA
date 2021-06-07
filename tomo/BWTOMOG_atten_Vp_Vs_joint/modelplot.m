addpath('matguts','plotting','function','seizmo','Tinycodes');
typelst = {'PZ_direct_num','PZ_direct_both','PZmix','Z_mag_5.5_10.0','P_mag_6.0_10.0'};
is2D = 0;
issave = 1;
for i = 1:length(typelst)
    clf;
    if ~is2D
        filename = ['model/model_',[typelst{i},'.mat']];
        titletxt = typelst{i};
    else
        filename = ['model/model_',[typelst{i},'_2D.mat']];
        titletxt = [typelst{i},'_2D'];
    end
    load(filename);
    %plot_Zmaps(plot_model,par,data,0);
    if issave
        plot_model = conv2plotable(model,par);
        plot_Hmaps(plot_model,par,1,par.saveopt)
        figure(67)
        print('-dpng',['/media/lun/easystore/tmp/poster/model/flatness/',titletxt,'_Hslice.png'])
        close(67)
        plot_Zmaps(plot_model,par,data,0)
        close(14)
        figure(15)
        print('-dpng',['/media/lun/easystore/tmp/poster/model/flatness/',titletxt,'_Zmap.png'])
        close(15)
    end
end