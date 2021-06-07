%{
typelst = {'PZmix','PZ_direct_num','P_mag_6.0_10.0','Z_mag_5.5_10.0'};
direct_lst = {'up','down'};
is2D_lst = {'','_2D'};
for k = 1:2
for j = 1:2
for it = 1:4

load(['sqztest/sqztest_',typelst{it},'_',direct_lst{k},is2D_lst{j},'.mat']);
for i = 1:length(squeezetest.zsqz)
    v(i) = squeezetest.freedom(i);
    var(i) = squeezetest.resid(i)^2/v(i);
end
disp([typelst{it},'_',direct_lst{k},is2D_lst{j}])
if strcmp(direct_lst{k},'up')
for i = length(squeezetest.zsqz):-1:2
    F = var(i)/var(i-1);
    p = 1 - fcdf(F,v(i),v(i-1));
    disp([num2str(par.zz(i)),'km vs ',num2str(par.zz(i-1)),'km: F=',num2str(F),', p=',num2str(p)])
end
else
for i = 1:length(squeezetest.zsqz)-1
    F = var(i)/var(i+1);
    p = 1 - fcdf(F,v(i),v(i+1));
    disp([num2str(par.zz(i)),'km vs ',num2str(par.zz(i+1)),'km: F=',num2str(F),', p=',num2str(p)])
end
end
fprintf('\n')
end
end
end
%}

issave = 1;
typelst = {'PZ_direct_num','PZmix','P_mag_6.0_10.0','Z_mag_5.5_10.0'};
direct_lst = {'up','down'};
is2D_lst = {'','_2D'};
for kt = 1:2 %1:2
for jt = 1:2 %1:2
for it = 1:4 %1:4

clf;
load(['sqztest/sqztest_',typelst{it},'_',direct_lst{kt},is2D_lst{jt},'.mat']);
p_lst = zeros(length(squeezetest.zsqz));
for i = 1:length(squeezetest.zsqz)
    v(i) = squeezetest.freedom(i);
    var(i) = squeezetest.resid(i)^2/v(i);
end
titletxt = [typelst{it},'_',direct_lst{kt},is2D_lst{jt}];
disp(titletxt)
if strcmp(direct_lst{kt},'up')
for i = length(squeezetest.zsqz):-1:1
    %for j = i-1:-1:1
    for j = length(squeezetest.zsqz):-1:1
    if i == j
        continue
    end
    F = var(i)/var(j);
    if F>1
        p = 1 - fcdf(F,v(i),v(j));
    else
        p = fcdf(F,v(i),v(j));
    end
    p_lst(i,j) = p;
    disp([num2str(par.zz(i)),'km vs ',num2str(par.zz(j)),'km: F=',num2str(F),', p=',num2str(p)])
    end
end
else
for i = 1:length(squeezetest.zsqz)
    for j = 1:length(squeezetest.zsqz)
    if i == j
        continue
    end
    %for j = i+1:length(squeezetest.zsqz)
    F = var(i)/var(j);
    if F>1
        p = 1 - fcdf(F,v(i),v(j));
    else
        p = fcdf(F,v(i),v(j));
    end
    p_lst(i,j) = p;
    disp([num2str(par.zz(i)),'km vs ',num2str(par.zz(j)),'km: F=',num2str(F),', p=',num2str(p)])
    end
end
end
fprintf('\n')

imagesc(p_lst'*100)
cmp = jet;
colormap(cmp)
clim = [5 50];
caxis(clim)
h = colorbar;
ylabel(h,'$Probability[\%]$','FontSize',12, 'Interpreter','Latex')
%colorbar;
for i = 1:length(squeezetest.zsqz)
    for j = 1:length(squeezetest.zsqz)
        if p_lst(i,j)>0 &&  p_lst(i,j)<=0.05
            rectangle('Position', [i-0.5, j-0.5, 1, 1], 'Edgecolor', 'none', 'FaceColor', 'black', 'LineWidth', 0.000001);
        elseif p_lst(i,j)==0
            rectangle('Position', [i-0.5, j-0.5, 1, 1], 'Edgecolor', 'white', 'FaceColor', 'white', 'LineWidth', 0.000001);
        end
    end
end

h1 = line(NaN,NaN,'LineWidth',5,'LineStyle','-','Color','black');
%h2 = line(NaN,NaN,'LineWidth',5,'LineStyle','-','Color','b');

for i = 1:length(squeezetest.zsqz)
    ticklabel{i} = num2str(squeezetest.zsqz(i));
end
%set(gca,'xaxislocation','top');
%set(gca,'xaxislocation','top');
axis equal
box()
xlabel('zsqz1[km]','interpreter', 'latex')
xticks(1:length(squeezetest.zsqz))
xticklabels(ticklabel)
xlim([0.5,length(squeezetest.zsqz)+0.5])
ylabel('zsqz2[km]','interpreter', 'latex')
yticks(1:length(squeezetest.zsqz))
yticklabels(ticklabel)
ylim([0.5,length(squeezetest.zsqz)+0.5])
grid()
set(gca,'ydir','normal')
title(['Ftest_',titletxt], 'interpreter', 'none')
legend([h1],{'$p \leq 5\%$'},'Interpreter','Latex')
if issave
    output = ['/media/lun/easystore/tmp/IM3.22/sqztest/Ftest/p_',titletxt,'.png'];
    print(gcf,'-dpng',output);
end
end
end
end