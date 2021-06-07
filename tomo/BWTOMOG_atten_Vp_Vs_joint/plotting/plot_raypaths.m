function plot_raypaths( par,data )
% plot_raypaths( par,data )
%   plot the raypaths, coloured by delay time
figure(82)
clf
hold on

load('rpath.mat')

times = data.ray.d;
timelims = mean(times)+[-3*std(times) 3*std(times)];
colo  = colour_get(times,max(timelims),min(timelims),jet);
[rpath.lat,rpath.lon] = project_xy(par,rpath.rx,rpath.ry,'inverse'); %#ok<NODEF>
% colour deets
caxis([min(times) max(times)])
% colormap = jet; %#ok<NASGU>

%plot
set(gcf,'DefaultAxesColorOrder',colo)
plot3(rpath.lon,rpath.lat,-rpath.rz)

% plot deets
title('Raypaths coloured by delay','FontSize',14)
xlabel('Longitude','FontSize',13);
ylabel('Latitude','FontSize',13);
xlim([min(par.mln) max(par.mln)])
ylim([min(par.mlt) max(par.mlt)])

% geographic data
addpath('/Users/zeilon/Dropbox/MATLAB/lib/borders_v3.1.2/borders/');
%borders('countries','k','linewidth',2)
% plot3(ncst(:,1),ncst(:,2),zeros(size(ncst(:,1))),'k','LineWidth',3) %#ok<NODEF>
% geoshow([S.Lat], [S.Lon], 'Color', 'black','linewidth',2)


% colour deets
caxis(timelims)
colormap(flipud(jet)); %#ok<NASGU>
colorbar;

set(gcf,'Position',[680   180   937   798])

view(-12,28);


end

