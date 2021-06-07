function [ datf ] = filt_quick( dat,flo,fhi,dt,npoles,Npass,dirpass,filtopt )
%[ datf ] = filt_quick( dat,fmin,fmax,dt,npoles,Npass,filtopt )
%   quick function to filter a trace (dat is Nx1 vector) or a set of traces
%   (dat is NxM matrix)
% 
%  if flo < ffund (fundamental freq = 1./T, where T is time series total
%  length) then this will apply a low pass using the high-f corner. 
%  if fhi == fNyq (samprate/2) then this will apply a high-pass using the
%  low-f corner. 
%  if both frequencies are outside of the ffund-fNyq range, then will just
%  return the data. 
%
%  Npass is the number 

if nargin<4
    dt = 1;
end
if nargin < 5
    npoles = 2;
end
if nargin < 6
    Npass = 2; % 
end
if nargin < 7
    dirpass = 1; % direction to run the filter if causal (forward = +1, backward = -1)
end
if nargin < 8
    filtopt = 'butter';
end

ffund = 1./(size(dat,1).*dt);

% calc the filter parms
switch filtopt
    case 'butter'
        % option 1: butter
        if flo > ffund && fhi.*dt.*2<1 % bandpass
            [bb,aa]=butter(npoles, [flo, fhi].*dt.*2);
        elseif flo > ffund && fhi.*dt.*2>=1 % high pass
            [bb,aa]=butter(npoles, flo.*dt.*2,'high');
        elseif  flo<=ffund && fhi.*dt.*2<1 % low pass
            [bb,aa]=butter(npoles, fhi.*dt.*2,'low');
        elseif  flo<=ffund && fhi.*dt.*2>=1
            datf = dat;
            return
        end
    case 'cheby'
        % option 2: cheby
        [bb,aa]=cheby1(npoles,0.5, [flo, fhi].*dt.*2);
    case 'highbutter'
        % option 3: higher order butter
        [z,p,k]=butter(npoles, [flo, fhi].*dt.*2.);
        [sos,g]=zp2sos(z,p,k); bb=sos; aa=g;
end


% do the filtering
datf = zeros(size(dat));
for ii = 1:size(dat,2)
    % pad
    dat_pad = [zeros(1000,1);dat(:,ii);zeros(1000,1)];
    
    if dirpass == -1, dat_pad = flipud(dat_pad); end
    
    % filter
    if Npass==1 % one-pass filter, causal, shape affected (dispersed)
        dat_padf=filter(bb, aa, dat_pad);
    elseif Npass==2 % two-pass, filter, non-causal (zero-phase), shape preserved
        dat_padf=filtfilt(bb, aa, dat_pad);
    end
    
    if dirpass == -1, dat_padf = flipud(dat_padf); end

    % unpad, store
    datf(:,ii) = dat_padf(1001:end-1000);
end


end

