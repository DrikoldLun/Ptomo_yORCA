datfile_lst = {'yORCA/data/PZmix.dat','yORCA/data/Z_mag_5.5_10.0.dat','yORCA/data/PZ_direct_num.dat','yORCA/data/P_mag_6.0_10.0.dat'};
direct_lst = {'up','down'};
is2D_lst = [0 1];

parfor i = 1:4
for j = 1:2
for k = 1:2
    squeezing_test(datfile_lst{i},direct_lst{j},is2D_lst(k));
end
end
end
