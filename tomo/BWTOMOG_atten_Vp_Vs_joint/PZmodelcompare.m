m1_lst = {'Z','PZmix','PZdirect','PZmix','PZdirect','PZmix'};
m2_lst = {'P','P','P','Z','Z','PZdirect'};
%m1_lst = {'Z','P','PZmix','PZdirect'};
depthrange = [0,300];
hitqthre = 0.3;
for j = 1:length(m1_lst)
clf;
m1 = m1_lst{j};
m2 = m2_lst{j};   
m = {m1,m2};
mval = {};
try
for i = 1:2
    %{
    if i == 1
    load(sprintf('flatmodel/model%s.mat',m1))
    elseif i == 2
    load(sprintf('smthmodel/model%s.mat',m1))
    end
    %}
    load(sprintf('flatmodel/model%s.mat',m{i}))
    if strcmp(m{i},'Z')
        mval{i} = modelZ.mval;
    elseif strcmp(m{i},'P') 
        mval{i} = modelP.mval;
    elseif strcmp(m{i},'PZmix')
        mval{i} = modelPZmix.mval;
    elseif strcmp(m{i},'PZdirect')
        mval{i} = modelPZdirect.mval;
    end
end
catch
disp('error, continue')
continue
end
         
mval1 = mval{1}(par.hitq>hitqthre&par.mz>depthrange(1)&par.mz<depthrange(2));
mval2 = mval{2}(par.hitq>hitqthre&par.mz>depthrange(1)&par.mz<depthrange(2));
plot(mval1,mval2,'o');
hold on;
plot([-1,1],[-1,1],'r')
axis equal
xlim([-0.03,0.03])
ylim([-0.03,0.03])
txt = sprintf("dep:%d-%dkm hitq\\_thre:%.1f",depthrange(1),depthrange(2),hitqthre);
title(txt)
xlabel(['model\_',m1])
ylabel(['model\_',m2])
%xlabel('flatness')
%ylabel('smoothness')
path = '/media/lun/easystore/tmp/poster/modelcompare/';
print(gcf,'-dpng',fullfile(path,sprintf('%s_%s.png',m1,m2)));
end