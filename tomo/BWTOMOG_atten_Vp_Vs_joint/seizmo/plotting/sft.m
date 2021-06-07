function [varargout]=sft(data,varargin)
%SFT    Sliding Fourier Transform of SEIZMO data
%
%    Usage:    sft(data)
%              sft(data,...,'width',percent,...)
%              sft(data,...,'units',type,...)
%              sft(data,...,'overlap',percent,...)
%              sft(data,...,'pow2pad',power,...)
%              sft(data,...,'freqrange',frange,...)
%              sft(data,...,'dbrange',dbrange,...)
%              sft(data,...,'dbcolormap',cmap,...)
%              sft(data,...,'axis',handle,...)
%              sft(data,...,'option',value,...)
%              data=sft(...)
%
%    Description:
%     SFT(DATA) plots records in SEIZMO struct DATA as the timeseries above
%     a spectrogram (for easy visual comparison).  Multiple records are
%     plotted in a series of subplots.  All records must be evenly sampled
%     time series records.  The spectrogram plot is generated by taking
%     windows of 2.5% the records' length with an overlap from one window
%     to the next of 75% for a total of 157 windows per record.  The
%     spectrogram is actually the power spectral density (PSD) in decibels
%     normalized so that 0dB is the maximum.  The plots are drawn in a new
%     figure window.
%
%     SFT(DATA,...,'WIDTH',PERCENT,...) specifies an alternative window
%     width.  Smaller windows will provide better time resolution at the
%     cost of frequency resolution.  The default width is 2.5%.
%
%     SFT(DATA,...,'UNITS',TYPE,...) changes the window width unit type.
%     The available choices are '%', 'N', & 'S'.  The default is '%'.  The
%     other choices allow for specifying the width in samples ('N') or in
%     seconds ('S').  Changing the type does not change the numeric value
%     of the window width (ie 2.5% becomes 2.5s).
%
%     SFT(DATA,...,'OVERLAP',PERCENT,...) indicates the overlap between
%     adjacent windows.  The default is 75%.  Higher overlap provides a
%     smoother spectrogram along the time axis.  Changing the window width
%     units does not change the overlap units!
%
%     SFT(DATA,...,'POW2PAD',POWER,...) controls the fft zero padding to
%     the next power of 2.  The default value is 0 which pads to the next
%     power of 2.  Higher powers provide higher frequency resolution but
%     can take considerably more computation.
%
%     SFT(DATA,...,'FREQRANGE',FRANGE,...) specifies the frequency range
%     to plot in the spectrogram.  The default FRANGE is [0 Fnyquist].
%     FRANGE must be [FREQLOW FREQHIGH].
%
%     SFT(DATA,...,'DBRANGE',DBRANGE,...) rescales the colormap to extend
%     to the decibel ranges in DBRANGE.  DBRANGE should be in [DBLO DBHI].
%     The default is [MINDB 0], which has no impact.
%
%     SFT(DATA,...,'DBCOLORMAP',CMAP,...) alters the colormap used in the
%     spectrogram plots.  The default colormap is FIRE.  The colormap may
%     be a Nx3 RGB triplet array or a string that may be evaluated to a Nx3
%     RGB triplet.
%
%     SFT(DATA,...,'AXIS',HANDLE,...) plots the entire set of spectrograms
%     in the space allocated to HANDLE.  Useful for compiling different
%     information into a single figure.
%
%     SFT(DATA,...,'OPTION',VALUE,...) sets certain plotting options to do
%     simple manipulation of the plots.  Available options are:
%      FGCOLOR    -- foreground color (axes, text, labels)
%      BGCOLOR    -- background color (does not set figure color)
%      COLORMAP   -- colormap for coloring records
%      XLABEL     -- record x axis label
%      YLABEL     -- y axis label
%      TITLE      -- title
%      XLIM       -- record x axis limits (tight by default)
%      YLIM       -- record y axis limits (tight by default)
%      LINEWIDTH  -- line width of records (default is 1)
%      LINESTYLE  -- line style of records (can be char/cellstr array)
%      NUMCOLS    -- number of subplot columns
%      UTC        -- plot in absolute time if TRUE (UTC w/ no leap support)
%      DATEFORMAT -- date format used if ABSOLUTE (default is auto)
%      XDIR       -- 'normal' or 'reverse'
%      YDIR       -- 'normal' or 'reverse'
%      FONTSIZE   -- size of fonts in the axes
%      FONTWEIGHT -- 'light', 'normal', 'demi' or 'bold'
%      MARKERS    -- true/false where true draws the markers
%
%     DATA=SFT(...) returns with the spectrogram PSD stored as the
%     dependent component data in SEIZMO struct DATA.  This means that the
%     records are of XYZ datatype.  B,E,DELTA,NXSIZE,XMINIMUM,XMAXIMUM
%     reflect window times.  NYSIZE, YMINIMUM, YMAXIMUM reflect the
%     frequency info.  NPTS is the number of pixels in the spectrogram.
%
%    Notes:
%
%    Header changes: B, E, DELTA, NPTS, NXSIZE, NYSIZE, IFTYPE,
%                    DEPMEN, DEPMIN, DEPMAX, XMINIMUM, XMAXIMUM,
%                    YMINIMUM, YMAXIMUM
%
%    Examples:
%     % Plot records with spectrograms showing 0-0.02Hz limited to -30dB:
%     sft(data,'fr',[0 0.02],'dbr',[-30 0]);
%
%     % Plot 4 records in a subplot:
%     h=subplot(3,3,5);
%     sft(data(1:4),'axis',h);
%
%    See also: DFT, CUT, SPECTROGRAM, PLOT1

%     Version History:
%        Apr. 26, 2010 - initial version
%        May   5, 2010 - fully working version
%        May   7, 2010 - fix overlap bug
%        July  8, 2010 - plot name fix & nargchk fix
%        Aug. 15, 2010 - linked x axes of record & spectragram, short idep
%        June 24, 2011 - change name from stft to sft to avoid conflict
%        Nov. 21, 2011 - bring up-to-date with plotting functions
%        Mar.  9, 2012 - drop getenum* usage
%        May  29, 2012 - pow2pad=0 by default, fix output breakage
%        Aug.  2, 2012 - set seizmocheck state back when no output, better
%                        checking (disallow empty for many options),
%                        transpose output to be consistent
%        Mar.  1, 2014 - maketime fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2014 at 18:25 GMT

% todo:
% - different spectrogram methods?
% - scaleogram?
% - non-logarithmic spectrum?

% check nargin
error(nargchk(1,inf,nargin));

% check struct
error(seizmocheck(data,'dep'));
nrecs=numel(data);

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% default/parse options
opt=parse_seizmo_plot_options(varargin{:});

% attempt spectrogram
try
    % check headers
    data=checkheader(data,...
        'FALSE_LEVEN','ERROR',...
        'XYZ_IFTYPE','ERROR');
    
    % verbosity
    verbose=seizmoverbose;
    
    % line coloring
    if(ischar(opt.CMAP) || iscellstr(opt.CMAP))
        % list of color names or a colormap function
        try
            % attempt color name to rgb conversion first
            opt.CMAP=name2rgb(opt.CMAP);
            opt.CMAP=repmat(opt.CMAP,ceil(nrecs/size(opt.CMAP,1)),1);
        catch
            % guess its a colormap function then
            opt.CMAP=str2func(opt.CMAP);
            opt.CMAP=opt.CMAP(nrecs);
        end
    else
        % numeric colormap array
        opt.CMAP=repmat(opt.CMAP,ceil(nrecs/size(opt.CMAP,1)),1);
    end
    
    % line style/width
    opt.LINESTYLE=cellstr(opt.LINESTYLE);
    opt.LINESTYLE=opt.LINESTYLE(:);
    opt.LINESTYLE=repmat(opt.LINESTYLE,...
        ceil(nrecs/size(opt.LINESTYLE,1)),1);
    opt.LINEWIDTH=opt.LINEWIDTH(:);
    opt.LINEWIDTH=repmat(opt.LINEWIDTH,...
        ceil(nrecs/size(opt.LINEWIDTH,1)),1);
    
    % check filetype (only timeseries or xy)
    iftype=getheader(data,'iftype id');
    spec=strcmpi(iftype,'irlim') | strcmpi(iftype,'iamph');
    
    % convert spectral to timeseries
    if(sum(spec)); data(spec)=idft(data(spec)); end
    
    % header info
    [b,npts,delta,ncmp,z6,kname,idep]=getheader(data,...
        'b','npts','delta','ncmp','z6','kname','idep desc');
    z6=datenum(cell2mat(z6));
    
    % get markers info
    [marknames,marktimes]=getmarkers(data);
    
    % convert markers to absolute time if used
    if(opt.ABSOLUTE)
        marktimes=marktimes/86400+z6(:,ones(1,size(marktimes,2)));
    end
    
    % valid string option values
    valid.UNITS={'%' 'per' 'percent' ...
                 'n' 'samp' 'samples' ...
                 's' 'sec' 'seconds'};
    
    % check sft options
    if(ischar(opt.UNITS)); opt.UNITS=cellstr(opt.UNITS); end
    if(~iscellstr(opt.UNITS) ...
            || any(~ismember(opt.UNITS,valid.UNITS)) ...
            || ~any(numel(opt.UNITS)==[1 nrecs]))
        error('seizmo:sft:badInput',...
            ['UNITS must be one of the following:' ...
            sprintf('%s ',valid.UNITS{:})]);
    end
    if(isscalar(opt.UNITS)); opt.UNITS(1:nrecs,1)=opt.UNITS; end
    if(~isreal(opt.WINDOW) || ~any(numel(opt.WINDOW)==[1 nrecs]) ...
            || any(opt.WINDOW<=0))
        error('seizmo:sft:badInput',...
            'WIDTH must be a real-valued scalar/vector >0!');
    end
    if(isscalar(opt.WINDOW)); opt.WINDOW(1:nrecs,1)=opt.WINDOW; end
    if(~isreal(opt.OVERLAP) || ~any(numel(opt.OVERLAP)==[1 nrecs]) ...
            || any(opt.OVERLAP<0 | opt.OVERLAP>100))
        error('seizmo:sft:badInput',...
            'OVERLAP must be between 0-100 (in %%)!');
    end
    if(isscalar(opt.OVERLAP)); opt.OVERLAP(1:nrecs,1)=opt.OVERLAP; end
    if(~isreal(opt.POW2PAD) || ~any(numel(opt.POW2PAD)==[1 nrecs]) ...
            || any(opt.POW2PAD~=fix(opt.POW2PAD)))
        error('seizmo:sft:badInput',...
            'POW2PAD must be an integer!');
    end
    if(isscalar(opt.POW2PAD)); opt.POW2PAD(1:nrecs,1)=opt.POW2PAD; end
    if(ischar(opt.DBCMAP))
        opt.DBCMAP=cellstr(opt.DBCMAP);
    end
    if(isreal(opt.DBCMAP) && ndims(opt.DBCMAP)==2 ...
            && size(opt.DBCMAP,2)==3 ...
            && all(opt.DBCMAP(:)>=0 & opt.DBCMAP(:)<=1))
        opt.DBCMAP={opt.DBCMAP};
    elseif(iscellstr(opt.DBCMAP))
        % nothing
    else
        error('seizmo:sft:badInput',...
            ['DBCOLORMAP must be a colormap function\n'...
            'string or a Nx3 RGB triplet array!']);
    end
    if(isscalar(opt.DBCMAP))
        opt.DBCMAP(1:nrecs,1)=opt.DBCMAP;
    elseif(numel(opt.DBCMAP)~=nrecs)
        error('seizmo:sft:badInput',...
            ['DBCOLORMAP must be a colormap function\n'...
            'string or a Nx3 RGB triplet array!']);
    end
    % secret option!!!
    if(~isreal(opt.HOTROD) ...
            || ~any(numel(opt.HOTROD)==[1 nrecs]) ...
            || any(opt.HOTROD<0 | opt.HOTROD>1))
        error('seizmo:sft:badInput',...
            'HOTROD must be between 0 & 1!');
    end
    if(size(opt.HOTROD,1)==1); opt.HOTROD(1:nrecs,1)=opt.HOTROD; end
    if(~isempty(opt.DBRANGE))
        if(~isreal(opt.DBRANGE) || size(opt.DBRANGE,2)~=2 ...
                || ~any(size(opt.DBRANGE,1)==[1 nrecs]) ...
                || any(opt.DBRANGE(:,1)>opt.DBRANGE(:,2)))
            error('seizmo:sft:badInput',...
                'DBRANGE must be Nx2 array of [DBLOW DBHIGH]!');
        end
        if(size(opt.DBRANGE,1)==1)
            opt.DBRANGE=opt.DBRANGE(ones(nrecs,1),:);
        end
    end
    if(~isempty(opt.FRANGE))
        if(~isreal(opt.FRANGE) || size(opt.FRANGE,2)~=2 ...
                || ~any(size(opt.FRANGE,1)==[1 nrecs]) ...
                || any(opt.FRANGE(:,1)>opt.FRANGE(:,2)) ...
                || any(opt.FRANGE<0))
            error('seizmo:sft:badInput',...
                'FREQRANGE must be Nx2 array of [FREQLO FREQHI]!');
        end
        if(size(opt.FRANGE,1)==1)
            opt.FRANGE=opt.FRANGE(ones(nrecs,1),:);
        end
    end
    
    % identify width units
    pw=ismember(opt.UNITS,{'%' 'per' 'percent'});
    nw=ismember(opt.UNITS,{'n' 'samp' 'samples'});
    sw=ismember(opt.UNITS,{'s' 'sec' 'seconds'});
    
    % convert width/overlap/pow2pad to proper units (samples)
    if(any(pw)); opt.WINDOW(pw)=ceil(opt.WINDOW(pw).*npts(pw)/100); end
    if(any(nw)); opt.WINDOW(nw)=ceil(opt.WINDOW(nw)); end
    if(any(sw)); opt.WINDOW(sw)=ceil(opt.WINDOW./delta); end
    opt.OVERLAP=min(opt.WINDOW-1,ceil(opt.OVERLAP.*opt.WINDOW/100));
    opt.POW2PAD=2.^(nextpow2n(opt.WINDOW)+opt.POW2PAD);
    
    % plotting setup
    if(~nargout)
        % error about plotting multi-cmp records
        if(any(ncmp>1))
            error('seizmo:sft:noMultiCMPplotting',...
                'SFT does not support plotting multicmp records!');
        end
        
        % make new figure if invalid or no axes handle passed
        if(isempty(opt.AXIS) || any(~ishandle(opt.AXIS)) ...
                || any(~strcmp('axes',get(opt.AXIS,'type'))))
            fh=figure('color',opt.BGCOLOR,'name','SFT -- SEIZMO');
            opt.AXIS=axes('parent',fh);
        end
        
        % handle single axis vs multiple axis
        if(numel(opt.AXIS)==nrecs)
            % get positioning
            drawnow;
            op=get(opt.AXIS,'outerposition');
            fh=get(opt.AXIS,'parent');
            delete(opt.AXIS);
            
            % uncell (if more than one record)
            if(iscell(op)); op=cell2mat(op); end
            if(iscell(fh)); fh=cell2mat(fh); end
            
            % "subplot positioning" aka not really
            pwidth=op(:,3);
            height=op(:,4);
            left=op(:,1);
            bottom=op(:,2);
        elseif(isscalar(opt.AXIS))
            % get positioning
            drawnow;
            op=get(opt.AXIS,'outerposition');
            fh=get(opt.AXIS,'parent');
            delete(opt.AXIS);
            
            % expand figure handle scalar so there is one for each record
            fh=fh(ones(nrecs,1));
            
            % number of columns
            if(isempty(opt.NUMCOLS))
                opt.NUMCOLS=round(sqrt(nrecs));
            end
            nrows=ceil(nrecs/opt.NUMCOLS);
            
            % outer position of each subplot
            % [left bottom width height]
            pwidth=op(3)/opt.NUMCOLS;
            height=op(4)/nrows;
            left=op(1):pwidth:op(1)+op(3)-pwidth;
            left=left(ones(nrows,1),:)';
            bottom=op(2):height:op(2)+op(4)-height;
            bottom=flipud(bottom(ones(opt.NUMCOLS,1),:)')';
        else
            error('seizmo:sft:badInput',...
                'Incorrect number of axes handles in AXIS!');
        end
        
        % plot constants to make text look "right"
        % ...well as much as we can without changing fontsize
        sh=0.6;
        sw=0.8;
        tpl=min([0.0875 (1-sh)*height sh*height sw*pwidth (1-sw)*pwidth]);
        jt=2/3*tpl;
        cbw=0.05;
    end
    
    % detail message
    if(verbose)
        disp('Getting Short Time Fourier Transform of Record(s)');
        print_time_left(0,nrecs);
    end
    
    % loop over records
    e=nan(nrecs,1); nxsize=e; nysize=e; depmen=e; depmin=e; depmax=e;
    for i=1:nrecs
        % skip dataless
        if(~npts(i)); continue; end
        
        % return spectrograms
        if(nargout)
            % get spectrogram
            P=cell(ncmp(i),1);
            for j=1:ncmp(i)
                [P{j},F,T,P{j}]=spectrogram(double(data(i).dep(:,j)),...
                    opt.WINDOW(i),opt.OVERLAP(i),...
                    opt.POW2PAD(i),1/delta(i));
                P{j}=P{j}.';
                P{j}=P{j}(:); % make column vector
            end
            
            % assign power spectra to dep
            data(i).dep=cell2mat(P);
            
            % get fields
            % Differences from SAC (!!!):
            % - b/e/delta are timing related
            % - xminimum/xmaximum/yminimum/ymaximum are time/freq values
            b(i)=b(i)+T(1);
            e(i)=b(i)+T(end);
            delta(i)=(e(i)-b(i))/(numel(T)-1);
            npts(i)=numel(P{1});
            nxsize(i)=numel(T);
            nysize(i)=numel(F);
            depmen(i)=nanmean(data(i).dep(:));
            depmin(i)=min(data(i).dep(:));
            depmax(i)=max(data(i).dep(:));
        else % plotting
            % how this is gonna look
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %        record name        %
            %    +---------------+-+    %
            % amp|  seismogram   |c|    %
            %    +---------------+b|    %
            %  f |               |a|    %
            %  r |  spectrogram  |r|    %
            %  e |               | |    %
            %  q +---------------+-+    %
            %        time (sec)   dB    %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % plot record
            figure(fh(i));
            h1=axes('parent',fh(i),...
                'outerposition',[left(i) bottom(i)+sh*height ...
                sw*pwidth (1-sh)*height],...
                'position',[left(i)+tpl bottom(i)+sh*height ...
                sw*pwidth-tpl (1-sh)*height-jt],...
                'ydir',opt.YDIR,'xdir',opt.XDIR,'color',opt.BGCOLOR,...
                'fontsize',opt.FONTSIZE,'fontweight',opt.FONTWEIGHT,...
                'xcolor',opt.FGCOLOR,'ycolor',opt.FGCOLOR);
            hold(h1,'on');
            if(opt.ABSOLUTE)
                rh=plot(h1,z6(i)+b(i)/86400+...
                    (0:delta(i)/86400:delta(i)/86400*(npts(i)-1)).',...
                    data(i).dep,...
                    'color',opt.CMAP(i,:),...
                    'linewidth',opt.LINEWIDTH(i),...
                    'linestyle',opt.LINESTYLE{i});
            else
                rh=plot(h1,b(i)+(0:delta(i):delta(i)*(npts(i)-1)).',...
                    data(i).dep,...
                    'color',opt.CMAP(i,:),...
                    'linewidth',opt.LINEWIDTH(i),...
                    'linestyle',opt.LINESTYLE{i});
            end
            hold(h1,'off');
            
            % let matlab sig proc box handle the psd estimation
            [P,F,T,P]=spectrogram(double(data(i).dep),...
                opt.WINDOW(i),opt.OVERLAP(i),opt.POW2PAD(i),1/delta(i));
            P=10*log10(abs(P));
            
            % set userdata to everything but data (cleared)
            data(i).dep=[];
            data(i).ind=[];
            data(i).index=[i 1];
            set(rh,'userdata',data(i));
            
            % tag records
            set(rh,'tag','record');
            
            % extras
            box(h1,'on');
            grid(h1,'on');
            axis(h1,'tight');
            
            % add markers to axis userdata
            userdata.markers.names=marknames(i,:);
            userdata.markers.times=marktimes(i,:);
            userdata.function='sft';
            set(h1,'userdata',userdata);
            
            % draw markers
            if(opt.MARKERS); drawmarkers(h1,varargin{:}); end
            
            % truncate to frequency range
            if(isempty(opt.FRANGE))
                frng=[min(F) max(F)];
            else
                frng=opt.FRANGE(i,:);
            end
            fidx=F>=frng(1) & F<=frng(2);
            F=F(fidx);
            P=P(fidx,:);
            
            % normalization
            maxp=max(P(:));
            P=P-maxp;
            
            % special normalization
            minp=min(P(:));
            P(P<minp*(1-opt.HOTROD(i)))=minp;
            
            % fix timing
            if(opt.ABSOLUTE)
                T=z6(i)+(T+b(i))/86400;
            else % relative
                T=T+b(i);
            end
            
            % deal with dbrange
            if(isempty(opt.DBRANGE))
                dbrng=[min(P(:)) max(P(:))];
            else
                dbrng=opt.DBRANGE(i,:);
            end
            
            % plot spectrogram ourselves
            h2=axes('parent',fh(i),...
                'outerposition',...
                [left(i) bottom(i) sw*pwidth sh*height],...
                'position',[left(i)+tpl bottom(i)+tpl ...
                sw*pwidth-tpl sh*height-tpl],...
                'ydir',opt.FDIR,'xdir',opt.XDIR,'color',opt.BGCOLOR,...
                'fontsize',opt.FONTSIZE,'fontweight',opt.FONTWEIGHT,...
                'xcolor',opt.FGCOLOR,'ycolor',opt.FGCOLOR);
            hold(h2,'on');
            imagesc(T,F,P,'parent',h2,dbrng);
            hold(h2,'off');
            
            % extras
            ylabel(h2,'Freq (Hz)');
            colormap(h2,opt.DBCMAP{i});
            axis(h2,'xy','tight');
            grid(h2,'on');
            
            % axis zooming
            if(~isempty(opt.XLIM)); xlim(h2,opt.XLIM); end
            if(~isempty(opt.YLIM)); ylim(h1,opt.YLIM); end
            
            % sync times of timeseries and spectrogram
            set(h1,'xlim',get(h2,'xlim'));
            
            % do the datetick thing
            if(opt.ABSOLUTE)
                if(isempty(opt.DATEFORMAT))
                    if(isempty(opt.XLIM))
                        datetick(h1,'x');
                        datetick(h2,'x');
                    else
                        datetick(h1,'x','keeplimits');
                        datetick(h2,'x','keeplimits');
                    end
                else
                    if(isempty(opt.XLIM))
                        datetick(h1,'x',opt.DATEFORMAT);
                        datetick(h2,'x',opt.DATEFORMAT);
                    else
                        datetick(h1,'x',opt.DATEFORMAT,'keeplimits');
                        datetick(h2,'x',opt.DATEFORMAT,'keeplimits');
                    end
                end
            else
                if(~isempty(opt.DATEFORMAT))
                    if(isempty(opt.XLIM))
                        datetick(h1,'x',opt.DATEFORMAT);
                        datetick(h2,'x',opt.DATEFORMAT);
                    else
                        datetick(h1,'x',opt.DATEFORMAT,'keeplimits');
                        datetick(h2,'x',opt.DATEFORMAT,'keeplimits');
                    end
                end
            end
            
            % turn off tick labels on timeseries
            set(h1,'xticklabel',[]);
            
            % now lets link the x axes together
            linkaxes([h1 h2],'x');
            
            % label
            if(~isempty(opt.TITLE) && isnumeric(opt.TITLE))
                switch opt.TITLE
                    case 1 % filename
                        if(~isempty(data(i).name))
                            p1title=texlabel(data(i).name,'literal');
                        else
                            p1title=['RECORD ' num2str(i)];
                        end
                    case 2 % kstnm
                        p1title=kname(i,2);
                    case 3 % kcmpnm
                        p1title=kname(i,4);
                    case 4 % shortcmp
                        p1title=kname{i,4}(3);
                    case 5 % stashort
                        p1title=strcat(kname(i,2),'.',kname{i,4}(3));
                    case 6 % stcmp
                        p1title=strcat(kname(i,2),'.',kname(i,4));
                    case 7 % kname
                        p1title=texlabel(strcat(kname(i,1),...
                            '.',kname(i,2),'.',kname(i,3),...
                            '.',kname(i,4)),'literal');
                    otherwise
                        p1title=['RECORD ' num2str(i)];
                end
            else
                p1title=opt.TITLE;
            end
            if(isnumeric(opt.XLABEL) && opt.XLABEL==1)
                if(opt.ABSOLUTE)
                    xlimits=get(h1,'xlim');
                    p1xlabel=joinwords(...
                        cellstr(datestr(unique(fix(xlimits)))),'   to   ');
                else
                    p1xlabel='Time (sec)';
                end
            else
                p1xlabel=opt.XLABEL;
            end
            if(isnumeric(opt.YLABEL) && opt.YLABEL==1)
                p1ylabel=idep(i);
            else
                p1ylabel=opt.YLABEL;
            end
            title(h1,p1title,'color',opt.FGCOLOR,...
                'fontsize',opt.FONTSIZE,'fontweight',opt.FONTWEIGHT);
            xlabel(h2,p1xlabel,'color',opt.FGCOLOR,...
                'fontsize',opt.FONTSIZE,'fontweight',opt.FONTWEIGHT);
            ylabel(h1,p1ylabel,'color',opt.FGCOLOR,...
                'fontsize',opt.FONTSIZE,'fontweight',opt.FONTWEIGHT);
            
            % plot colorbar
            % entire height
            % 25% of width
            c=colorbar('peer',h2,'position',...
                [left(i)+sw*pwidth bottom(i)+tpl ...
                cbw*pwidth height-tpl-jt]);
            xlabel(c,'dB');
            set(c,'xcolor',opt.FGCOLOR,'ycolor',opt.FGCOLOR);
            set(get(c,'XLabel'),'color',opt.FGCOLOR);
        end
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end
    
    % update headers if there is output args
    if(nargout)
        varargout{1}=changeheader(data,'b',b,'e',e,'delta',delta,...
            'npts',npts,'iftype','ixyz','nxsize',nxsize,'nysize',nysize,...
            'depmen',depmen,'depmin',depmin,'depmax',depmax,...
            'xminimum',b,'xmaximum',e,...
            'yminimum',0,'ymaximum',1./(2*delta));
    end
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
