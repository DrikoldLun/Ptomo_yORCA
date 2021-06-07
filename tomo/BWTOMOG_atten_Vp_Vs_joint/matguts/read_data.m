function [ data,par ] = read_data( DataFile,StaFile,CrustFile,par )

nph = length(par.phases);
%{
% check data type
for ip = 1:nph
if par.t_ts == 1
    if ~any(regexp(DataFile{ip},'dT'))
        error('Datatype indicated as dT ==> maybe wrong data file')
    end
elseif par.t_ts == 2
    if ~any(regexp(DataFile{ip},'dtstar'))
        error('Datatype indicated as dtstar ==> maybe wrong data file')
    end
end
end
%}

%=================== input rays =======================================
clear dinfo;
if ~iscell(DataFile), DataFile = {DataFile}; end
for ip = 1:nph
    disp(DataFile{ip});
    fid = fopen(DataFile{ip}, 'r');
    %A = textscan(fid, '%7.4f  %7.3f  %7.2f  %s  %s  %4.0f  %8.4f  %9.4f  %7.3f  %6.2f  %6.4f  %6.4f');
    A = textscan(fid, '%*s %f %f %f %s %s %f %f %f %f %f %*f %f %*f %*f %f %*s');
    fclose(fid);
    % how many fields of A
    nflds = size(A,2); 
    % how many obs
    nrays = size(A{1},1);
    % add info about phase (P=1, S=2)
    A{nflds+1} = find(strcmp({'P','S'},par.phases{ip}(1)))*ones(nrays,1);
    nflds = size(A,2); 
    
    if ~exist('dinfo','var')
        dinfo = A;
    else
        % augment orids
        addorid = 10.^ceil(log10(max(dinfo{6})));
        A{6} = addorid + A{6};
        for icol = 1:nflds
            dinfo{icol} = cat(1,dinfo{icol},A{icol});
        end
    end
    clear A
        
        
end
    

ray.pd    = dinfo{1};
ray.p     = r2d(ray.pd)/6371;
ray.gcarc = dinfo{2};
%ray.baz   = mod(dinfo{3}-par.map_proj.origin(3),360);
ray.baz   = dinfo{3};
ray.sta   = dinfo{4};
ray.nwk   = dinfo{5};
ray.orid  = dinfo{6};
ray.elat  = dinfo{7};
ray.elon  = dinfo{8};
ray.edep  = dinfo{9};
ray.d     = dinfo{10};
ray.sd    = dinfo{11};
ray.cf    = dinfo{12};
ray.ph    = dinfo{13};
clear dinfo

ray.nrays = length(ray.p);

data.ray=ray;

%======================= save event info ==============================

[evt.orid, I] = unique(ray.orid);
evt.lat   = ray.elat(I);
evt.lon   = ray.elon(I);
evt.depth = ray.edep(I);
evt.nevts = length(evt.orid);

data.evt = evt;

%====================== input station =================================
[sta_names,ia] = unique( ray.sta );
sta_nwks = ray.nwk(ia);

stn.alllat = [];
stn.alllon = [];
if ~iscell(StaFile), StaFile = {StaFile}; end
for ip = 1:nph
    fid = fopen(StaFile{ip}, 'r');
    B = textscan(fid,'%s %s %f %f %f ','delimiter',','); 
    fclose(fid);
    nflds = size(B,2);
    if ~exist('sinfo','var')
        sinfo = B;
    else % concat
        % augment orids
        for icol = 1:nflds
            sinfo{icol} = cat(1,sinfo{icol},B{icol});
        end
    end
    stn.alllat = cat(1,stn.alllat,B{3})
    stn.alllon = cat(1,stn.alllon,B{4})
    clear B 
end
% find unique STA/NW codes
[iunq] = stanw_unique(sinfo{1},sinfo{2});


% stas{is}, slats(is), slons(is), selevs(is), statype{is}, IP,trm,...
%         xrdg,ocage,plate,yr,moh);

stns.sta = strtok(sinfo{1}(iunq),' '); % remove trailing spaces
stns.nwk = strtok(sinfo{2}(iunq),' '); % remove trailing spaces
stns.lat = sinfo{3}(iunq);  
stns.lon = sinfo{4}(iunq);  
stns.elv = sinfo{5}(iunq)/1000;  % convert to km
%stns.elv = -4.62588*ones(length(stns.lon),1);
% stns.typ = strtok(B{5},' ');  % remove trailing spaces
% stns.ip  = B{6}; 
% stns.trm = B{7};  
% stns.Xrd = B{8};
% stns.age = B{9};
% stns.plt = strtok(B{10},' ');  % remove trailing spaces
% stns.yr  = B{11};
% stns.moh = B{12};
% stns.sed = B{13}/1000; % go from m to km
clear sinfo
% stns.sed(stns.sed==-.999) = 0;


stn.lat = zeros(length(sta_names),1);
stn.lon = zeros(length(sta_names),1);
stn.elv = zeros(length(sta_names),1);
stn.moh = zeros(length(sta_names),1);
stn.vpvs = zeros(length(sta_names),1);
stn.nrf = zeros(length(sta_names),1);
% stn.typ = cell(length(sta_names),1);
% stn.ip  = cell(length(sta_names),1);
% stn.trm = zeros(length(sta_names),1);
% stn.Xrd = zeros(length(sta_names),1);
% stn.age = zeros(length(sta_names),1);
% stn.plt = cell(length(sta_names),1);
% stn.yr  = zeros(length(sta_names),1);
% stn.moh = zeros(length(sta_names),1);
% stn.sed = zeros(length(sta_names),1);

%====================== input crust =================================
stns.moh = nan(length(stns.lat),1);
stns.vpvs = nan(length(stns.lat),1);
stns.nrf = nan(length(stns.lat),1);
if ~isempty(CrustFile) && exist(CrustFile,'file')
    fid = fopen(CrustFile,'r');
    C = textscan(fid,'%s%s%f%f%f%f%f','delimiter',','); 
    fclose(fid);
    
    % parse into stn structure
    for is = 1:length(stns.lat)
        indx = strcmp(C{1},stns.sta{is}) & strcmp(C{2},stns.nwk{is});
        if ~any(indx), continue; end
        stns.moh(is) = C{5}(indx);
        stns.vpvs(is) = C{6}(indx);
        stns.nrf(is) = C{7}(indx);
    end
    

end


% ------- associate data to stations --------------------------------------------
for ii = 1:length(sta_names)
    try
    ix=find(strcmp(sta_names{ii}, stns.sta ) & strcmp(sta_nwks{ii}, stns.nwk ));
    
    stn.lat(ii) = stns.lat(ix);
    stn.lon(ii) = stns.lon(ix);
    stn.elv(ii) = stns.elv(ix);
    stn.moh(ii) = stns.moh(ix);
    stn.vpvs(ii) = stns.vpvs(ix);
    stn.nrf(ii) = stns.nrf(ix);
    catch 
        1
    end
%     stn.typ(ii) = stns.typ(ix);
%     stn.ip(ii)  = stns.ip(ix);
%     stn.trm(ii) = stns.trm(ix);
%     stn.Xrd(ii) = stns.Xrd(ix);
%     stn.age(ii) = stns.age(ix);
%     stn.plt(ii) = stns.plt(ix);
%     stn.yr(ii)  = stns.yr(ix);
%     stn.moh(ii) = stns.moh(ix);
%     stn.sed(ii) = stns.sed(ix);
    stn.num(ii) = ii;
end

stn.sta = sta_names;
stn.nstas = length(sta_names);

ray.sta_num = zeros(length(ray.d),1);
for jj = 1:length(stn.sta)
    ix=find( strcmp(stn.sta{jj},ray.sta) );
    ray.sta_num(ix) = jj;
    ray.stalat(ix,1) = stn.lat(jj);
    ray.stalon(ix,1) = stn.lon(jj);
end


% data structure
data.ray = ray;
data.stn = stn;
par.stn = stn;

% add the estimate of the total ray lengths
data = trdist_est(data, par);

end

%% subfunctions

function [ data ] = trdist_est( data, par )
    %crude estimate works well enough
    % these rayp-distance values are from taup, using IASP91

    ipdp   = [ 0     4.548   4.638   5.403  6.148  6.877  7.602  8.304  8.615  8.835  9.099  10.898 13.627 13.700];
    ipds   = [ 0     8.658   9.199   10.519 11.724 12.868 13.963 14.958 15.368 15.669 15.961 20.047 24.314 24.560];
    idist     = [ 13000 10563.5 10007.5 8895.6 7783.6 6672   5560   4448   3892   3336   2780   2224   1668   1112];

    pd = data.ray.pd;
    
    if par.PS==1
        pray = true(data.ray.nrays,1);
        sray = false(data.ray.nrays,1);
    elseif par.PS==2
        pray = false(data.ray.nrays,1);
        sray = true(data.ray.nrays,1);
    elseif par.PS==3  
        pray = data.ray.ph==1;
        sray = data.ray.ph==2;
    end
    
    data.ray.trdist = nan(data.ray.nrays,1);
    % do p
    data.ray.trdist(pray) = interp1(ipdp,idist,pd(pray));
    % do s
    data.ray.trdist(sray) = interp1(ipds,idist,pd(sray));

    
    
    
end
