issave = 1;
typelst = {'PZ_direct_num','PZmix','P_mag_6.0_10.0','Z_mag_5.5_10.0'};
is2D_lst = {'','_2D'};
for jt = 1:1 %1:2
for it = 1:1 %1:4

clf;
load(['sqztest_2side/sqz2side_',typelst{it},is2D_lst{jt},'.mat']);
post_E = ones(length(par.zz),length(par.zz))*0;
titletxt = [typelst{it},is2D_lst{jt}];
disp(titletxt)
fprintf('\n')
nc = length(squeezetest.is);
index = unique(squeezetest.is);
depth = unique(squeezetest.zs);
for ic = 1:nc
    post_E(squeezetest.is(ic),squeezetest.id(ic)) = squeezetest.wvr(ic)/squeezetest.freedom(ic);%/squeezetest.freedom(ic);
end

imagesc(post_E')
cmp = jet;
colormap(cmp)
clim = [max(post_E(:))*0.7 max(post_E(:))];
caxis(clim)
h = colorbar;
ylabel(h,'Normalized wvr[\%]','FontSize',12, 'Interpreter','Latex')
[i_op,j_op] = find(post_E==max(post_E(:)));
hold on;
plot(i_op,j_op,'o','MarkerEdgeColor','none','MarkerSize',10,'MarkerFaceColor','white')
%colorbar;
for i = index
    for j = index
        if i>j
            rectangle('Position', [i-0.5, j-0.5, 1, 1], 'Edgecolor', 'white', 'FaceColor', 'white', 'LineWidth', 0.000001);
        end
    end
end

%h1 = line(NaN,NaN,'LineWidth',5,'LineStyle','-','Color','black');
%h2 = line(NaN,NaN,'LineWidth',5,'LineStyle','-','Color','b');

for i = index
    ticklabel{i} = num2str(depth(i));
end
%set(gca,'xaxislocation','top');
%set(gca,'xaxislocation','top');
axis equal
box()
xlabel('$Z_{shallow}[km]$','interpreter', 'latex')
xticks(index)
xticklabels(ticklabel)
xlim([min(index)-0.5,max(index)+0.5])
ylabel('$Z_{deep}[km]$','interpreter', 'latex')
yticks(index)
yticklabels(ticklabel)
ylim([min(index)-0.5,max(index)+0.5])
grid()
set(gca,'ydir','normal')
title(['sqz2side_',titletxt], 'interpreter', 'none')
%legend([h1],{'p<=0.05'})
if issave
    output = ['/media/lun/easystore/tmp/IM3.22/sqztest/Econtour_sqz2side_',titletxt,'_wvrreso.png'];
    print(gcf,'-dpng',output);
end
end
end