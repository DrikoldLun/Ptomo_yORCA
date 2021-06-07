function [ array,ofile ] = array_geom(centpt,az,pod,opt,kmlfile)
% [ array ] = array_geom(option)
% 
% Function to spit out the array geometry, in both array-based x,y
% coordinates and the lat/lon coordinates. 
% 
% INPUTS:
%   centpt - centre of the array [lat, lon]
%   az     - azimuth of array (of x=0 line), degrees from N
%   kmlfile - file name for kml output
%   option - option for deployment configuration. Default is 1.

if nargin<5 
    kmlfile = '';
end

if nargin<5 || isempty(opt)
    opt=1;
end

nstas = 30;
for is = 1:nstas, stas{is,1} = ['OBS_',num2str(is)]; end

switch opt % different deployment options
    case 1 % simple grid with dense line down the x=0 line
        sxx = []; syy = [];
        sxx = [sxx; -250*ones(4,1)]; syy = [syy; (-187.5:125:187.5)']; % far left line
        sxx = [sxx; -125*ones(5,1)]; syy = [syy; (-200:100:200)']; % mid left line
        sxx = [sxx;  0*ones(12,1)]; syy = [syy; (-192.5:35:192.5)']; % dense centre line
        sxx = [sxx;  125*ones(5,1)]; syy = [syy; (-200:100:200)']; % mid right line
        sxx = [sxx;  250*ones(4,1)]; syy = [syy; (-187.5:125:187.5)']; % far right line
    case 2 % concentric circles + centre line + corners
        sxx = []; syy = [];
        sxx = [sxx; 200*sin(2*pi*(0:7)'/8)]; syy = [syy; 200*cos(2*pi*(0:7)'/8)]; % 200 km circle of 8
        sxx = [sxx; 250*sin(2*pi*(0.5:7.5)'/8)]; syy = [syy; 250*cos(2*pi*(0.5:7.5)'/8)]; % 250 km off-set circle of 8
        sxx = [sxx; 100*sin(2*pi*(0:5)'/6)]; syy = [syy; 100*cos(2*pi*(0:5)'/6)]; % 100 km circle of 6
        sxx = [sxx;  0*ones(4,1)]; syy = [syy; [-150;-100/3;100/3;150]]; % dense centre line
        sxx = [sxx;  [-250;-250;250;250]]; syy = [syy; [-250;250;-250;250]]; % corners
    case 3 % Don's tweak to simple grid with dense line down the x=0 line
        sxx = []; syy = [];
        sxx = [sxx; -250*ones(5,1)]; syy = [syy; (-200:100:200)']; % far left line
        sxx = [sxx; -125*ones(4,1)]; syy = [syy; (-187.5:125:187.5)']; % mid left line
        sxx = [sxx;  0*ones(12,1)]; syy = [syy; (-192.5:35:192.5)']; % dense centre line
        sxx = [sxx;  125*ones(4,1)]; syy = [syy; (-187.5:125:187.5)']; % mid right line
        sxx = [sxx;  250*ones(5,1)]; syy = [syy; (-200:100:200)']; % far right line
        sxx(sxx~=0) = sxx(sxx~=0) + random('norm',0,5,size(sxx(sxx~=0)));
        syy(sxx~=0) = syy(sxx~=0) + random('norm',0,5,size(syy(sxx~=0)));
    case 4 % Jim's tweak to extend dense line down the x=0 line
        sxx = []; syy = [];
        sxx = [sxx; -250*ones(3,1)]; syy = [syy; (-140:140:140)']; % far left line
        sxx = [sxx; -125*ones(5,1)]; syy = [syy; (-200:100:200)']; % mid left line
        sxx = [sxx;  0*ones(14,1)]; syy = [syy; (-227.5:35:227.5)']; % dense centre line
        sxx = [sxx;  125*ones(5,1)]; syy = [syy; (-200:100:200)']; % mid right line
        sxx = [sxx;  250*ones(3,1)]; syy = [syy; (-140:140:140)']; % far right line
        sxx(sxx~=0) = sxx(sxx~=0) + random('norm',0,10,size(sxx(sxx~=0)));
        syy(sxx~=0) = syy(sxx~=0) + random('norm',0,10,size(syy(sxx~=0)));
        sxx(sxx==0) = sxx(sxx==0) + random('norm',0,1,size(sxx(sxx==0)));
        syy(sxx==0) = syy(sxx==0) + random('norm',0,1,size(syy(sxx==0)));
    case 5 % simple 4554 grid with variable spaced dense line down the x=0 line
        sxx = []; syy = [];
        sxx = [sxx; -250*ones(4,1)]; syy = [syy; (-187.5:125:187.5)']; % far left line
        sxx = [sxx; -125*ones(5,1)]; syy = [syy; (-200:100:200)']; % mid left line
        sxx = [sxx;  0*ones(12,1)]; syy = [syy; [-225:48:-129,-92:36.9:93,129:48:225]']; % dense centre line
        sxx = [sxx;  125*ones(5,1)]; syy = [syy; (-200:100:200)']; % mid right line
        sxx = [sxx;  250*ones(4,1)]; syy = [syy; (-187.5:125:187.5)']; % far right line
    case 6 % tweaked 4554 grid with variable dense line down the x=0 line
        sxx = []; syy = [];
        sxx = [sxx; [-250 -225 -225 -250]']; syy = [syy; (-187.5:125:187.5)']; % far left line
        sxx = [sxx; [-125 -105 -125 -105 -125]']; syy = [syy; (-200:100:200)']; % mid left line
        sxx = [sxx;  0*ones(12,1)]; syy = [syy; [-225:48:-129,-92:36.9:93,129:48:225]']; % dense centre line
        sxx = [sxx;  [125 105 125 105 125]']; syy = [syy; (-200:100:200)']; % mid right line
        sxx = [sxx;  [250 225 225 250]']; syy = [syy; (-187.5:125:187.5)']; % far right line
    case 7 % tweaked 4554 grid with variable dense line down the x=0 line
        sxx = []; syy = [];
        sxx = [sxx; [-250 -225 -225 -250]']; syy = [syy; (-187.5:125:187.5)']; % far left line
        sxx = [sxx; [-125 -105 -125 -105 -125]']; syy = [syy; (-180:90:180)']; % mid left line
        sxx = [sxx; [0;normrnd(0,1,10,1);0]]; syy = [syy; [-225:48:-129,-92:36.9:93,129:48:225]']; % dense centre line
        sxx = [sxx;  [125 105 125 105 125]']; syy = [syy; (-180:90:180)']; % mid right line
        sxx = [sxx;  [250 225 225 250]']; syy = [syy; (-187.5:125:187.5)']; % far right line
end


%% Project to map coords
proj.origin  = centpt;   % model origin
mstruct = defaultm('mercator');
mstruct.origin = [proj.origin -az];
mstruct = defaultm( mstruct ); 
proj.map_proj = mstruct;

[slats, slons] = project_xy(proj, sxx, syy, 'inverse');

%% station elevations
%{
if strcmp(pod,'Npod')
    topofile = '~/Dropbox/Work/Proposals/2016_PacificArray/plotting/etopo1_Npod.grd';
elseif strcmp(pod,'Spod')
    topofile = '~/Dropbox/Work/Proposals/2016_PacificArray/plotting/etopo1_Spod.grd';
end
[topoX,topoY,topoZ] = grdread2(topofile); % load topo grid
[topoX,topoY] = meshgrid(double(topoX),double(topoY)); topoZ = double(topoZ);
selevs = interp2(topoX,topoY,topoZ,slons,slats)/1000; % in units of km
if any(isnan(selevs))
    error('Maybe loaded wrong topo grid')
end
%}
%% Spit out array details


array = struct('stas',{stas},'nstas',nstas,...
                 'sxx',sxx,'syy',syy,...
                 'slats',slats,'slons',slons,...
                 'proj',proj,...
                 'pod',pod,'opt',opt);


%% write kml file
if ~isempty(kmlfile)
    kmldir = '/Users/zeilon/Work/Proposals/2016_PacificArray/plotting/';
    kmlwritepoint([kmldir,kmlfile],slats,slons,'name','')
end

%% write station data file
ofile = ['stafile_',pod,num2str(opt),'.dat'];
fid = fopen(ofile,'w');
fprintf(fid,'sta   ,  lat    ,  lon     ,  elev  ,  x    , y\n') ;
for is = 1:nstas
    fprintf(fid,'%-6s, %8.4f, %9.4f, %6.1f, %6.1f\n',...
        char(stas(is)), slats(is), slons(is),sxx(is),syy(is));
end
fclose(fid) ;   

%% plots
ifplot = 0;
if ifplot

figure(77), clf
plot(sxx,syy,'o','Markersize',10,'MarkerFaceColor','b')
axis([-290 290 -220 220])
grid on
axis equal


figure(78), clf
plot(slons,slats,'o','Markersize',10,'MarkerFaceColor','b')
axis equal

end

end

