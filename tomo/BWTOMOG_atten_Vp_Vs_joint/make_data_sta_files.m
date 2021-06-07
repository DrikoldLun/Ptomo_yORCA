%% parse results into large nstas*nevts structures
clear all
% if running alone, establish these, otherwise use previous values
phase = 'SKS';
component = 'R';
method = 'specR';

%% conditions
mag_min = 6.25; % skip events if below this mag
acor_min = 0.8; % skip events if xcor max below this acor
snr_min = 10; % skip result if data below this snr

%% parms
dTsd = 0.05; % standard minimum standard deviation. Divide this by acor
dtssd = 0.1; % standard minimum standard deviation. Divide this by acor


%% directories 
% project details
dbname = 'EARdb';
dbdir = '/Users/zeilon/Work/EastAfrica/EARdb/'; % include final slash

%% Preliminaries
addpath('matguts')
run([dbdir,dbname,'_startup.m']);

tomodatdir = ['/Users/zeilon/Dropbox/MATLAB/Atten_tomog_rays/',dbname,'/data/'];

%% =================================================================== %%
%% ==========================  GET TO WORK  ========================== %%
%% =================================================================== %%

% event details
load([infodir,'/events'],'evinfo'); 
% station details
load([infodir,'/stations'],'stainfo');

%% out files
dT_file = sprintf('data_dT_%s%s.dat',phase,component);
dtstar_file = sprintf('data_dtstar_%s_%s%s.dat',method,phase,component);
sta_file = sprintf('stations_%s%s.dat',phase,component);



%% db data
% [ nstas,stas,slats,slons,selevs,ondate,offdate,staname,statype ] = db_stadata( dbdir,dbnam );
% [ norids,orids,elats,elons,edeps,evtimes,mags ]  = db_oriddata( dbdir,dbnam );
% sages = jdf_crust_age(slats,slons,'linear');

%% =================================================================== %%
%% ============================  STAS  =============================== %%
%% =================================================================== %%

% fields: sta  slat  slon  selev (moho)
fid = fopen([tomodatdir,sta_file],'w');
for is = 1:stainfo.nstas
    fprintf(fid,'%-6s, %-6s, %8.4f, %9.4f, %7.1f \n',...
        stainfo.stas{is}, stainfo.nwk{is}, stainfo.slats(is), stainfo.slons(is), stainfo.selevs(is));
end
fclose(fid) ;   

%% =================================================================== %%
%% ===========================  DATA  ============================== %%
%% =================================================================== %%
dT_data = struct('rayp',[],'gcarc',[],'seaz',[],'sta',{},'nwk',{},...
                  'orid',[],'elat',[],'elon',[],'edep',[],...
                  'dat',[],'sd',[],'cfreq',[]);
dtstar_data = struct('rayp',[],'gcarc',[],'seaz',[],'sta',{},'nwk',{},...
                  'orid',[],'elat',[],'elon',[],'edep',[],...
                  'dat',[],'sd',[],'cfreq',[]);
              

count_dT = 0;
count_dts = 0;
              
for ie = 1:evinfo.norids % 44:norids % loop on orids
    fprintf('Orid %.0f\n',ie)
    if  evinfo.evmags(ie)<mag_min, continue, end
    orid = evinfo.orids(ie);
    
    %% grab data!
    % name files and directories
    evdir = [num2str(orid,'%03d'),'_',evinfo.datestamp(ie,:),'/'];
    datinfofile = [datadir,evdir,'_datinfo_',phase];
    arfile      = [datadir,evdir,'_EQAR_',phase,'_',component];
    
    % check files exist
    if ~exist([datinfofile,'.mat'],'file'), fprintf('No station mat files for this event\n');continue, end
    if ~exist([arfile,'.mat'],'file'), fprintf('No arfile for this event + phase\n');continue, end    
    % load info file
    load(datinfofile) % loads datinfo stucture
    % options to skip
    if isempty(datinfo), fprintf('No station mat files for this event\n'); continue, end
    if ~isfield(datinfo,'xcor'), fprintf('   NEED TO XCOR\n',phase), continue, end

    % load data file
    load(arfile)      % loads eqar structure
    if all([eqar.pred_arrT]==0), fprintf('No %s arrivals for this event\n',phase), continue, end
    
%% parse dT data
    if ~isfield(eqar,'dT'), fprintf('   NEED TO XCOR\n',phase), continue, end
    if length([eqar.dT])~=length(eqar),for is = 1:length(eqar),if isempty(eqar(is).dT),eqar(is).dT=nan;end,end,end
    dTgind = find(~isnan([eqar.dT]));
    dTgind = dTgind([eqar(dTgind).acor] >= acor_min);
    mean_dT = mean([eqar(dTgind).dT]);
    for is = 1:length(dTgind)
        ind = dTgind(is);
        count_dT = count_dT+1;
        
        dT_data(count_dT,1) ...
         = struct('rayp',eqar(ind).rayp,...
                  'gcarc',eqar(ind).gcarc,...
                  'seaz',eqar(ind).seaz,...
                  'sta',eqar(ind).sta,'nwk',which_sta_nwk(stainfo,eqar(ind).sta,eqar(ind).slat,eqar(ind).slon),...
                  'orid',orid,'elat',evinfo.elats(ie),'elon',evinfo.elons(ie),'edep',evinfo.edeps(ie),...
                  'dat',eqar(ind).dT - mean_dT,...
                  'sd',dTsd./eqar(ind).acor,...
                  'cfreq',mean([eqar(ind).par_dT.fhi,eqar(ind).par_dT.flo]));
    end
    
    if strcmp(method,'comb')
        
        %% parse dtstar_COMB data
        if ~isfield(eqar,'dtstar_comb'), fprintf('   NEED TO CALC DTSTAR w COMB\n',phase), continue, end
        if length([eqar.dtstar_comb])~=length(eqar), error('Not enough dtstars somehow'), end
        dtsgind = find(~isnan([eqar.dtstar_comb]));
        mean_dtstar = mean([eqar(dtsgind).dtstar_comb]);
        for is = 1:length(dtsgind)
            ind = dtsgind(is);
            count_dts = count_dts+1;

            dtstar_data(count_dts,1) ...
             = struct('rayp',eqar(ind).rayp,...
                      'gcarc',eqar(ind).gcarc,...
                      'seaz',eqar(ind).seaz,...
                      'sta',eqar(ind).sta,'nwk',which_sta_nwk(stainfo,eqar(ind).sta,eqar(ind).slat,eqar(ind).slon),...
                      'orid',orid,'elat',evinfo.elats(ie),'elon',evinfo.elons(ie),'edep',evinfo.edeps(ie),...
                      'dat',eqar(ind).dtstar_comb - mean_dtstar,...
                      'sd',dtssd,...
                      'cfreq',mean([1./eqar(ind).par_dtstar_comb.comb.Tmin,1./eqar(ind).par_dtstar_comb.comb.Tmax]));
        end
    
    elseif strcmp(method,'specR')
        
        %% parse dtstar_specR data
        if ~isfield(eqar,'dtstar'), fprintf('   NEED TO CALC DTSTAR w specR\n',phase), continue, end
        if length([eqar.dtstar])~=length(eqar), error('Not enough dtstars somehow'), end
        dtsgind = find(~isnan([eqar.dtstar]));
        mean_dtstar = mean([eqar(dtsgind).dtstar]);
        for is = 1:length(dtsgind)
            ind = dtsgind(is);
            count_dts = count_dts+1;

            dtstar_data(count_dts,1) ...
             = struct('rayp',eqar(ind).rayp,...
                      'gcarc',eqar(ind).gcarc,...
                      'seaz',eqar(ind).seaz,...
                      'sta',eqar(ind).sta,'nwk',which_sta_nwk(stainfo,eqar(ind).sta,eqar(ind).slat,eqar(ind).slon),...
                      'orid',orid,'elat',evinfo.elats(ie),'elon',evinfo.elons(ie),'edep',evinfo.edeps(ie),...
                      'dat',eqar(ind).dtstar - mean_dtstar,...
                      'sd',dtssd,...
                      'cfreq',mean([eqar(ind).fcrosslo,eqar(ind).fcrosshi]));
        end
    end

end % loop on orids

fprintf('%.0f dtstar measurements for %s on %s comp\n',count_dts,phase,component)
fprintf('%.0f dT measurements for %s on %s comp\n',count_dT,phase,component)

%% =================================================================== %%
%% =============================  dT  ================================ %%
%% =================================================================== %%

% fields: rayp  gcarc  seaz  sta  nwk orid  elat  elon  edep  dT  sd  c_freq
fid = fopen([tomodatdir,dT_file],'w');
for ii = 1:count_dT
fprintf(fid,'%7.4f  %7.3f  %7.2f  %6s  %6s  %4.0f  %8.4f  %9.4f  %7.3f  %6.2f  %6.4f  %6.4f\n',...
    dT_data(ii).rayp,dT_data(ii).gcarc,dT_data(ii).seaz,...
    dT_data(ii).sta,dT_data(ii).nwk,...
    dT_data(ii).orid,dT_data(ii).elat,dT_data(ii).elon,dT_data(ii).edep,...
    dT_data(ii).dat,dT_data(ii).sd,dT_data(ii).cfreq);
end
fclose(fid);
   


%% =================================================================== %%
%% ===========================  dtstar  ============================== %%
%% =================================================================== %%

% fields: rayp  gcarc  seaz  sta  nwk orid  elat  elon  edep  dtstar  sd  c_freq
fid = fopen([tomodatdir,dtstar_file],'w');
for ii = 1:count_dts
fprintf(fid,'%7.4f  %7.3f  %7.2f  %6s  %6s  %4.0f  %8.4f  %9.4f  %7.3f  %6.2f  %6.4f  %6.4f\n',...
    dtstar_data(ii).rayp,dtstar_data(ii).gcarc,dtstar_data(ii).seaz,...
    dtstar_data(ii).sta,dtstar_data(ii).nwk,...
    dtstar_data(ii).orid,dtstar_data(ii).elat,dtstar_data(ii).elon,dtstar_data(ii).edep,...
    dtstar_data(ii).dat,dtstar_data(ii).sd,dtstar_data(ii).cfreq);
end
fclose(fid);

return
%#ok<*UNRCH>
%% =================================================================== %%
%% ==========================  MERGE S  ============================== %%
%% =================================================================== %%

%% dtstar
ST_dtstar_file = [tomodatdir,'data_dtstar_',method,'_ST.dat']; 
SKSR_dtstar_file = [tomodatdir,'data_dtstar_',method,'_SKSR.dat']; 
Sall_dtstar_file = [tomodatdir,'data_dtstar_',method,'_Sall.dat']; 
% read ST
fid = fopen(ST_dtstar_file,'r');
A_st = textscan(fid,'%7.4f  %7.3f  %7.2f  %6s  %6s  %4.0f  %8.4f  %9.4f  %7.3f  %6.2f  %6.4f  %6.4f\n');
fclose(fid);
% find max ST orid
maxoridst = max(A_st{6});

% read SKSR
fid = fopen(SKSR_dtstar_file,'r');
A_sksr = textscan(fid,'%7.4f  %7.3f  %7.2f  %6s  %6s  %4.0f  %8.4f  %9.4f  %7.3f  %6.2f  %6.4f  %6.4f\n');
fclose(fid);
% augment min SKSR orid
A_sksr{6} = A_sksr{6} + maxoridst;

% write Sall
for icol = 1:size(A_st,2)
    A_sall{icol} = [A_st{icol};A_sksr{icol}];
end
count_dtstar = length(A_sall{1});
% fields: rayp  gcarc  seaz  sta  nwk orid  elat  elon  edep  dtstar  sd  c_freq
fid = fopen(Sall_dtstar_file,'w');
for ii = 1:count_dtstar
fprintf(fid,'%7.4f  %7.3f  %7.2f  %6s  %6s  %4.0f  %8.4f  %9.4f  %7.3f  %6.2f  %6.4f  %6.4f\n',...
    A_sall{1}(ii),A_sall{2}(ii),A_sall{3}(ii),A_sall{4}{ii},A_sall{5}{ii},A_sall{6}(ii),...
    A_sall{7}(ii),A_sall{8}(ii),A_sall{9}(ii),A_sall{10}(ii),A_sall{11}(ii),A_sall{12}(ii));
end
fclose(fid);

%% dT
ST_dT_file = [tomodatdir,'data_dT_ST.dat']; 
SKSR_dT_file = [tomodatdir,'data_dT_SKSR.dat']; 
Sall_dT_file = [tomodatdir,'data_dT_Sall.dat']; 
% read ST
fid = fopen(ST_dT_file,'r');
A_st = textscan(fid,'%7.4f  %7.3f  %7.2f  %6s  %6s  %4.0f  %8.4f  %9.4f  %7.3f  %6.2f  %6.4f  %6.4f\n');
fclose(fid);
% find max ST orid
maxoridst = max(A_st{6}); 

% read SKSR
fid = fopen(SKSR_dT_file,'r');
A_sksr = textscan(fid,'%7.4f  %7.3f  %7.2f  %6s  %6s  %4.0f  %8.4f  %9.4f  %7.3f  %6.2f  %6.4f  %6.4f\n');
fclose(fid);
% augment min SKSR orid
A_sksr{6} = A_sksr{6} + maxoridst;

% write Sall
for icol = 1:size(A_st,2)
    A_sall{icol} = [A_st{icol};A_sksr{icol}];
end
count_dtstar = length(A_sall{1});
% fields: rayp  gcarc  seaz  sta  nwk orid  elat  elon  edep  dT  sd  c_freq
fid = fopen(Sall_dT_file,'w');
for ii = 1:count_dtstar
fprintf(fid,'%7.4f  %7.3f  %7.2f  %6s  %6s  %4.0f  %8.4f  %9.4f  %7.3f  %6.2f  %6.4f  %6.4f\n',...
    A_sall{1}(ii),A_sall{2}(ii),A_sall{3}(ii),A_sall{4}{ii},A_sall{5}{ii},A_sall{6}(ii),...
    A_sall{7}(ii),A_sall{8}(ii),A_sall{9}(ii),A_sall{10}(ii),A_sall{11}(ii),A_sall{12}(ii));
end
fclose(fid);

%% station file
copyfile([tomodatdir,sta_file],[tomodatdir,'stations_Sall.dat'])

    
    
    
    