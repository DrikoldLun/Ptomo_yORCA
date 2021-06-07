%OBS_DO_ROTATION
%  Do all the rotations for all the OBS data to get into proper Z-N-E form
%  Input:   needs portapng db  with .wfdisc .site .origin  table
%           needs origin table to only have events which satisfy correct
%           mag,distance conditions
%  Output:  portapng_rot.wfdisc & data in data/rotSTA directory 
% ZJE 03/2012
% orids=[1:21]; ?????? change
%   OBS names and corrections 
$$$$$$$ MAKE SURE YOU BACK UP THE WFDISC!
stas= {'B'};    
%All:   {'B','D', 'E', 'F', 'G', 'H', 'J'};
%Land:  {'PEMM', 'KEIA', 'JONE', 'GOGO'};  
scors= [12.1]; 
%All:   [ 116.8  0.3  183.1  301  47.2  324.5  143.9]; 
%Land:  [15.6  12.1  21.4  17.2]; 

dbnam='pngbodydb';
dbout='pngbodydb_rot';

windstart = 0.1; %seconds after EQ to open window
windend = 1799.9;   %seconds after EQ to close window
ichans=char({'BHE','BHN','BHZ'}); % input channels in order e ,n ,z
ochans=char({'BHE','BHN','BHZ'}); % output channels in order e ,n ,z
fchans=char({'e','n','z'});       % channel suffices for sac files

% Processing starts here

%loop over stations and records
for is=1:length(stas)
sta=char(stas(is));
odir=sprintf('data/rot%s/',sta);
xdir=sprintf('data/pre_rot%s',sta);
scor=scors(ksta);
sscor=sind(scor);
cscor=cosd(scor);
    
db=dbopen(dbnam,'r');
dbwf    = dblookup_table(db,'wfdisc');
dbsite  = dblookup_table(db,'site');
dbsitechan = dblookup_table(db,'sitechan');
dbjws=dbjoin(dbwf,dbsite);
dbs1 = dbsubset(dbwf,sprintf('sta =="%s"',sta));
nrecs = dbnrecs(dbs1);
wfids = dbgetv(dbs1,'wfid');
usedwfids = zeros(size(wfids));

% loop through records, trying to find all three components for each time,
% and cycling up in increasing wfid. 
for irec = 1:nrecs
wfid = wfids(irec);
if any(usedwfids == wfid)==1, continue; end % move on if already done this wfid
fprintf('Rotating sta %s, record %u/%u\n',sta,irec,nrecs);

%% find all other records starting coincident with this one
dbs1.record = dbfind(dbs1,sprintf('wfid == %u',wfid));
rtime = dbgetv(dbs1,'time');
dbs1.record=-501; %reset to whole table
dbs2 = dbsubset(dbs1,sprintf('time - %e < 1 && %e - time > 0 ',rtime+0.5,rtime+0.5)); %subset to same record start time
t0s = dbgetv(dbs2,'time');
t1s = dbgetv(dbs2,'endtime');
t0 = max(t0s); %starttime that will include all records
t1 = min(t1s); %endtime that will include all records

%% Check all chans are there
dbs2e = dbsubset(dbs2,sprintf('chan == %s',ichans(1)));
dbs2n = dbsubset(dbs2,sprintf('chan == %s',ichans(2)));
dbs2z = dbsubset(dbs2,sprintf('chan == %s',ichans(3)));
if dbnrecs(dbs2e)<1

%% Only excerpt if both the horiz chans are there
%% Old East
[tte, dat_e, ~, nsamps, samprate, wfide] = dbgetwfz(dbs2e, sta, t0, t1, 'epoch', ichans(1));
%% Old North
[ttn, dat_n, ~, nsamps, samprate, wfidn] = dbgetwfz(dbs2n, sta, t0, t1, 'epoch', ichans(2));
%% Old Vertical
[ttz, dat_z, ~, nsamps, samprate, wfidz] = dbgetwfz(dbs2z, sta, t0, t1, 'epoch', ichans(3));

$ need to do something if either horiz is not there = still remove from db
if isempty(dat_e)==1 || isempty(dat_n)==1
    

$ do some dbfinds to now assign each channel to a record
$ at the end, have to add the used wfids to the usedwfids
$ give new chans the old chans' wfids
$ move the old data into a "pre_rot" folder
$ place new data into "rot" db and folder
break


% view pre-rotated
for ic=1:3
    figure(13);
    subplot(3,2,2*ic-1); plot(tt,dat(:,ic),'b');
    title(sprintf('Component %s pre',char(chans(ic))));
end

% ROTATE 
        data_N =  data_n.*cscor - data_e.*sscor;
        data_E =  data_n.*sscor + data_e.*cscor;
        data_Z =  data_z;
        dat=[data_E,data_N,data_Z];
        
% view post-rotated
for ic=1:3
    figure(13);
    subplot(3,2,2*ic); plot(tt,dat(:,ic),'r');
    title(sprintf('Component %s post',ochans(ic,:)));
end

%% Put into wfdisc 
yn='nall';
for ic=1:3
% Set up trace and its parameters
tr=trnew;
tr=dblookup_table(tr,'trace');
wfchan=ochans(ic,:);
dfile=sprintf('%s.%s.sac.%s',epoch2str(evtime,'%Y.%j.%H.%M.%S'),sta,fchans(ic));

% delete old file query?
ofile=dir(sprintf('%s/%s',odir,dfile));
if isempty(ofile)==0 && strcmp(yn,'yall')==0
    yn=input('Do you want to overwrite file? (y/n/yall) ','s');
end
if strcmp(yn,'y')==1 || strcmp(yn,'yall')==1 
    delete(sprintf('%s/%s',odir,dfile));
elseif strcmp(yn,'n')==1
    continue
end

% Put the waveform into the trace-object:
tr.record=dbaddv(tr,'net', 'ROT','sta',sta,'chan',wfchan,...
'nsamp',nsamps,'samprate',samprate,'time',tt0,'endtime', tt1);
trinsert_data(tr,dat(:,ic));
% Save the trace data in a new database, with the underlying file in sac format:
dbo=dbopen(dbout,'r+');
dbow=dblookup_table(dbo,'wfdisc');
trsave_wf(tr,dbow,'sc',sprintf('%s/%s',odir,dfile));
trdestroy( tr ) % Clean up
dbclose(dbo)
end % loop on channels
end % loop on orids
end % loop on stations

% dbclose(db);
