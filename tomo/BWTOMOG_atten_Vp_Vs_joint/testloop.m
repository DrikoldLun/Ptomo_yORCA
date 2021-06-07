%datfile = 'yORCA/data/data_dT_PP.dat';
datfile = 'yORCA/data/PZmix.dat'; %2,5
%datfile = 'yORCA/data/Z_mag_5.5_10.0.dat'; %1,5
%datfile = 'yORCA/data/PZ_direct_num.dat'; %2,5
%datfile = 'yORCA/data/P_mag_6.0_10.0.dat'; %2,5
%datfile = 'yORCA/data/PZ_direct_both.dat'; %2,5
synth_test_ORCA(datfile,2,5);
%{
for damp = 0:1:5
    for smooth = 0:1:5
        synth_test_ORCA(datfile,damp,smooth);
    end
end
%}

%{
datfile_lst = dir('yORCA/data/syn*_dT_PZ.dat');

for i = 1:length(datfile_lst)
if strcmp(datfile_lst(i).name,'syn3998_dT_PZ.dat')
    continue
end
datfile = fullfile(datfile_lst(i).folder,datfile_lst(i).name);
for damp = 0:1:5
    for smooth = 0:1:5
        synth_test_ORCA(datfile,damp,smooth);
    end
end
end
%}