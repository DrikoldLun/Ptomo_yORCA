function [ data ] = makewaves_fn(db,windopt,t0,t1,orid,filtf1,filtf2,ndec,assocopt)
% [ data ] = makewaves_fn(db,windopt,t0,t1,orid,filtf1,filtf2,ndec,assocopt)
% This function makes a large structure of all the data across the array for
% a given arrival
%
% bit like dbgetwfX but for a whole array.
%
% For a given database, and a ref location, this script will filter events
% by distance and then output all possible orids
% For a chosen orid, the script will then loop through stations and load
% all the waveforms associated with the orid and put them into the
% cell "alldata" which has nsta + 1 structures:
% the first nsta structures contain the info and data for each station
% the last structure contains the evinfo
%
% INPUTS
% db = an open antelope database, with wfdisc, origin and site tables
% windopt   = 'time' or name of a phase for arrival mode
% t0        = beginning of time window or time before arrival
% t1        = end of time window or time after arrival
% orid      = have to know orid number in origin table
% filtfreqs = filtering frequencies
% ndec      = factor to decimate by, if at all
% assocopt  = '1' means only do stas with assoc for this orid. Default is 0
% 
% OUTPUTS
% data = structure with nstas layers
%           -first layer has 
%               orid
%               phase
%             	elat
%               elon
%               edep
%               evtime
%           -each layer has
%               sta
%           	slat
%           	slon
%               selev
%               seaz
%           	foraz
%               gcarc
%           	rayParam
%           	tt0 -- vector of absolute times (epochal)
%           	tte -- vector of times (s) since orid
%               tdiff -- diff between pred and pick (+ive = early)
%               samprate
%           	nsamps
%           	dof
%           	dt
%           	datN
%               datE
%               datR -- +ive away from source
%           	datT -- +ive to LEFT looking away from source
%           	datZ
%% put in goodorids function

%% starting info
% rlat=-10;
% rlon=150;
buffer=1; % for dbgetwfz wfdisc join
samprate=50;
corephases={'PKS','SKS','SKKS'};
if nargin < 8; ndec=1; end
if nargin < 9; assocopt=0; end
%% windowing: 
% excerpts by time
windstart = t0; %seconds after EQ to open window
windend = t1;   %seconds after EQ to close window
%OR around arrival
phase=windopt; % phase to window around
pretime=t0;
posttime=t1;

if strcmp(phase,'time')==1, windl=t1-t0-2*buffer;
else windl=t0+t1;
end


%% get all event info
dbsite=dblookup_table(db,'site');
dbor=dblookup_table(db,'origin');
dbwf=dblookup_table(db,'wfdisc');
dbassoc=dblookup_table(db,'assoc');

%% set orid and loop through stations getting data
dboo=dbor;
dboo.record=dbfind(dbor,sprintf('orid==%d',orid));
   
elat=dbgetv(dboo,'lat');
elon=dbgetv(dboo,'lon');
edep=dbgetv(dboo,'depth');
evtime=dbgetv(dboo,'time');
dbjos=dbjoin(dbsite,dbor);
dbjows=dbjoin(dbjos, dbwf,{'sta', 'sarrival()'},{'sta','time::endtime'});
dbjowso=dbsubset(dbjows,sprintf('orid==%d',orid));
if assocopt==1 && strcmp(phase,'time')~=1;
    dbjowso = dbjoin(dbjowso,dbassoc,{'sta','orid'},{'sta','orid'});
    dbjowso = dbsubset(dbjowso,sprintf('phase=="%s"',phase));
end
stas=unique(dbgetv(dbjowso,'sta'));
nstas=length(stas);
if nstas==0, error('No traces...'); end

%% set up structure for all data - station data/info goes in the first nstas rows and
% the last row is for a structure with the orid information
data=struct('orid',orid,...
            'phase',windopt,...
            'elat',elat,...
            'elon',elon,...
            'edep',edep,...
            'evtime',evtime,...
          'sta',cell(1,1),...
          'slat',zeros(1,1),...
          'slon',zeros(1,1),...
          'selev',zeros(1,1),...
          'seaz',zeros(1,1),...
          'foraz',zeros(1,1),...
          'gcarc',zeros(1,1),...
          'rayParam',zeros(1,1),...
          'tt0',zeros(1,1),...
          'tte',zeros(1,samprate*windl),...
          'tdiff',zeros(1,1),...
          'samprate',zeros(1,1),...
          'nsamps',zeros(1,1),...
          'dof',zeros(1,1),...
          'dt',zeros(1,1),...
          'datN',zeros(1,samprate*windl/ndec),...
          'datE',zeros(1,samprate*windl/ndec),...
          'datR',zeros(1,samprate*windl/ndec),...
          'datT',zeros(1,samprate*windl/ndec),...
          'datZ',zeros(1,samprate*windl/ndec));
      
% Making these this big seems unecessary - CHECK
for is=1:nstas
sta=char(stas(is));
dbsta=dbsubset(dbsite,sprintf('sta=="%s"',sta));
slat=dbgetv(dbsta,'lat');
slon=dbgetv(dbsta,'lon');
selev=dbgetv(dbsta,'elev');
seaz=azimuth(slat,slon,elat,elon);
foraz=azimuth(elat,elon,slat,slon);
gcarc=distance(slat,slon,elat,elon);
% move on if too close for good core phases
if any(strcmp(phase,corephases))==1 && gcarc < 85
    data(is).sta='mark';
    fprintf('Too close for discernible core phases\n')
    continue
end

fprintf('Processing sta %s orid %d phase %s filt [%.3f %.3f]\n',sta,orid,phase,filtf1,filtf2);

if strcmp(windopt,'time')==1, 
    t0=evtime+windstart+buffer;
    t1=evtime+windend-buffer;
else 
    tt=taupTime([],edep,phase,'sta',[slat,slon],'evt',[elat elon]);
    if isempty(tt)==1, data(is).sta='mark'; % mark for deletion
    continue, end % if no arrival
    arrtime=tt(1).time;
    rayParam=tt(1).rayParam;
    t0=evtime+arrtime-pretime;
    t1=evtime+arrtime+posttime;
end

%% get data
try
[tt, dat,chans,nsamps,samprate, wfids] = dbgetwfz(db,sta,t0,t1,'epoch');
catch me
chans=0;
end
chanend=char({'e','n','z'});
nchan=size(chans,1);
% only continue if all  channels exist
if nchan~=3; 
    data(is).sta='mark'; % mark for deletion
    continue
end % if enough chans

%% Data info
dat=detrend(dat);
ampl=max(max(abs(dat)));
dat=dat ./ ampl;
%%%% taper - window...
dt = 1/samprate;
tte = tt - evtime; % time since event
% nchan=length(chans);
saz=sin(pi-seaz*pi/180); % pi shift to rotate so radial is in foraz
caz=cos(pi-seaz*pi/180);

% Things I want in the structure...
data(is).sta=sta;
data(is).slat=slat;
data(is).slon=slon;
data(is).selev=selev;
data(is).seaz=seaz;
data(is).foraz=foraz;
data(is).gcarc=gcarc;
if strcmp(windopt,'time')~=1,
data(is).arrtime=arrtime;
data(is).rayParam=rayParam;
if assocopt==1 
dbs1 = dbsubset(dbassoc,sprintf('sta == "%s" && orid == %u && phase == "%s"',sta,orid,phase));
data(is).tdiff = dbgetv(dbs1,'timeres');
end
end
cd /Users/Zach/Documents/MATLAB/PNG_swsplit/

%% Data processing on three channels
for k=1:nchan
eval(sprintf('data_%s=dat(:,k);',chanend(k)));
end

% FILTER %%% filt filt?
[b,a]=butter(3,[filtf1,filtf2].*2.*dt);
data_e=filter(b,a,data_e);
data_n=filter(b,a,data_n);
data_z=filter(b,a,data_z);
% data_e=filtfilt(b,a,data_e);
% data_n=filtfilt(b,a,data_n);
% data_z=filtfilt(b,a,data_z);
amp2=max([max(abs(data_e)),max(abs(data_n)),max(abs(data_z))]);
data_e=data_e./amp2;
data_n=data_n./amp2;
data_z=data_z./amp2;
tt0=tt;

% WINDOW - taper
we=window(@tukeywin,length(data_e));
wn=window(@tukeywin,length(data_n));
wz=window(@tukeywin,length(data_z));
data_e=data_e.*we;
data_n=data_n.*wn;
data_z=data_z.*wz;

%  DECIMATE here after filtering
if (ndec > 1)
  nold=length(data_z);
  rsmp=1:ndec:nold;
  tt0=tt0(rsmp);
  dtold=dt;
  dt=mean(diff(tt0));
  data_e=data_e(rsmp);
  data_n=data_n(rsmp);
  data_z=data_z(rsmp);
  nsamps=nsamps/ndec;
end

% ROTATE (T positive to the left!) (recall caz+saz have been rotated so
% angle is the forward azimuth of the arrival
data_r =  caz.*data_n + saz.*data_e; %#ok<NASGU>
data_t =  saz.*data_n - caz.*data_e; %#ok<NASGU>

chan=char({'N','E','R','T','Z'});
chanend=char({'n','e','r','t','z'});
nchan=size(chan,1);
dof=zeros(nchan,2); % work out degrees of freedom
for k=1:nchan % put all traces into data structure
eval(sprintf('data(is).dat%s=data_%s;',chan(k),chanend(k)));
eval(sprintf('[dof(k,1),dof(k,2)]=scdofcalc(data_%s);',chanend(k)));% work out degrees of freedom
end
%assign to output structure
% eval(sprintf('%s.twind=[windstart-buffer windend+buffer];',sta));

data(is).samprate=samprate;
data(is).nsamps=nsamps;
data(is).tt0=tt0;
data(is).tte=tt0-evtime;
data(is).dt=dt;
data(is).dof=dof;
% %reset channels
% chan=char({'e','n','z'});
% nchan=size(chan,1);
% assign to alldata
% eval(sprintf('alldata{j}=%s;',sta));
% eval(sprintf('clear(''%s'');',sta));
end % loop on stations
%% delete null rows
marked=0;
for is=1:nstas
    if strcmp(data(is).sta,'mark')==1, marked(length(marked)+1)=is; end
end
data(marked(2:end))=[]; clear('marked');

end
