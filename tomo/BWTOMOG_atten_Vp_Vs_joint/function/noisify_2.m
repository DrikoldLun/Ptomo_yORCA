function   [ odat,SNRest,dofps ] = noisify_2( idat,ndat,samprate,SNR,smoothing,plotopt,theta,phi,dT )
% function [ odat,SNRest,dofps ] = noisify_2( idat,ndat,samprate,SNR,smoothing,plotopt,theta,phi,dT )
% 
% Function to take an input set of data (idat) and add noise with the power
% distribution of an input noise series (ndat). Ideally the two input
% series should be of the same length and sample rate
%
% INPUTS:
% idat       - time series of data in columns (e.g. 3 channels would be an
%                nsamps x nchans matrix) (try to have E,N,Z)
% ndat       - time series of several noisy data clips in columns 
%                i.e. nsamps x nchans x nclips Noise window should have no arrival. 
% samprate   - sample rate of input time series (samples per second)
% SNR        - signal to noise ratio desired
% smoothing  - 1 to smooth noise psd, anything else will not.
% plotopt    - 1 to plot, anything else will not.
% -- if you have them --
% theta      - forward azimuth of incoming SKS wave (relative to north)
% phitrue    - true fast direction (relative to north)
% dttrue     -  true time delay between fast and slow
% 
% OUTPUTS:
% odat       - output time series, same size as idat and with noise added
% SNRest     - measured SNR of odat
% dofps      - degrees of freedom per second of the output series
%
% Define SNR as max(initial signal)/2*std(noise) the max of the initial
% signal is not the same as the max value of the post-split signal (idat)
% as some of the energy will be on the T component, assuming some splitting
% Thus the SNRin, above will have to be the desired SNR * max(idat)/signal
% where signal is the pre-split max
% 
%
global psde

if nargin < 5
    smoothing=0;
end
if nargin < 6
    plotopt=0;
end

%% Key values
noise=max(max(abs(idat)))./SNR;
dt=1./samprate;
nsamps=size(idat,1); 
nchans=size(idat,2);
nclips=size(ndat,3);
samplength=nsamps/samprate;
% must be as many chans in noise and noise length must be integer multiple
% of the data length, so frequencies match
if size(ndat,2)~=nchans || size(ndat,1)/nsamps ~= round(size(ndat,1)/nsamps) 
    error('Make ndat the right size! Must have nchans and n*nsamps')
end

%% Detrend and scale input time series
for ii = 1:nclips 
ndat(:,:,ii) = detrend(ndat(:,:,ii),'constant'); %detrend
% idat=detrend(idat,'constant'); %detrend
%%%%%%%%% adjust amplitudes %%%%%%%%%%%%%%
for ic=1:nchans; ndat(:,ic,ii)=ndat(:,ic,ii)./sqrt(mean(ndat(:,ic,ii).^2)); end 
% for ic=1:nchans; idat(:,ic)=idat(:,ic)./sqrt(mean(idat(:,ic).^2)); end 
end

%% How much of input is signal?
    for ic=1:nchans
    indsig = find(abs(idat(:,ic)) > 0); Lsig(ic) = length(indsig);
    end
%% GET PSD of NOISE
% Window and get sample of noise
% get psd for each channel
if isempty(psde)~=1
    disp('N.B. Using existing psd')
else
nw=4;
psde=zeros(nsamps/2 + 1,nchans+1,nclips);
freq=[1:floor(nsamps/2)+1]'./(nsamps*dt);%in mHz - to get in Hz, divide by 1000
for ii = 1:nclips
for ic=1:nchans
    [Pxx,~] = pmtm(ndat(1:nsamps,ic,ii),nw,freq,samprate);
    H = dspdata.psd(Pxx,'Fs',samprate);
    if plotopt==1
        figure(10)
        subplot(nchans,nclips+1,ii + (ic-1)*(nclips+1))
        plot(H)
        ylim([-150,50])
    end
    psde(:,1,ii)=H.Frequencies;
    psde(:,ic+1,ii)=H.Data;
    if smoothing==1
        psde(:,ic+1,ii) = smooth(psde(:,ic+1,ii),round(samprate/3));
        if plotopt==1
            figure(10);	subplot(nchans,nclips+1,ii + (ic-1)*(nclips+1))
            hold on;	plot(H); ylim([-150,50]);	hold off
        end
    end
end %loop on chans
end % loop on clips
end % if using existing psd

% Average across clips to get representative psd of noise
psde = mean(psde,3);
if plotopt==1
    for ic=1:nchans
        figure(10); subplot(nchans,nclips+1,(nclips+1) + (ic-1)*(nclips+1)); hold on
        plot(psde(:,1),10*log10(psde(:,ic+1)),'r'); ylim([-150,50])
        title('Power Spectral Density')
        xlabel('Frequency (Hz)')
        ylabel('Power/frequency (dB/Hz)'); grid on
    end
end

%% Put in noise 
% make fft of synth traces and add to fft of a bootstrap noise series
%%%% REPEAT OVER NOISE REALISATIONS UNTIL SNRest IS CLOSE TO SNR IN %%%%
dSNR=1000; 
while dSNR > SNR/10 
%     disp('go')
clear('i');
FA=zeros(nsamps,nchans);
FB=zeros(nsamps,nchans);
%% Ft of noise
for ic=1:nchans
FA(1,ic)=psde(1,ic+1);
FA(nsamps/2+1,ic)=psde(end,ic+1);
FA(2:nsamps/2,ic)=psde(2:end-1,ic+1).*exp(1i*2*pi*random('unif',0,1,nsamps/2-1,1));
FA(nsamps/2+2:end,ic)=flipud(conj(FA(2:nsamps/2,ic)));
end
%% Ft of signal
for ic=1:nchans
FB(:,ic)=fft(idat(:,ic));
end

%% Have to think carefully about powers and amplitude scaling
Pi = sum(idat.^2,1); % get power of input so can scale to unit power
Pif = sum(abs(FB).^2,1); % get power of input Ft so can scale to unit power
Pnf = sum(abs(FA).^2,1);% get power of initial FA so can scale to unit power
Phorz = mean(Pnf(1:2)); % average power on horizontal channels
for ic = 1:nchans
% % power in the fourier domain is sum of squares - so scale amps by sqrt(1/P)
% FA(:,ic)=FA(:,ic)/sqrt(Pnf(ic)); % FA power now normalised to 1 on each channel
FA(:,ic)=FA(:,ic)/sqrt(Phorz); % FA power now normalised to ~1 on each horizontal
% % ifft will cause factor of 1/N change in power - so scale amps by sqrt(N)
FA(:,ic)=FA(:,ic)*sqrt(nsamps);
% % Finally must scale to 10*power of idat, as noise series is 10xlonger than signal pulse
% FA(:,ic)=FA(:,ic) * sqrt(10*sum(Pi));
% % Finally must scale to (N/4)*power of peak idat to get SNR = peak/2*std
FA(:,ic)=FA(:,ic) * sqrt(nsamps/4);
% FB(:,ic)=fft(idat(:,ic)/A(ic)); % for now do not scale by amplitude of original here
end
% Now the power on each of the horizontals should be equal to 
%((N * max(R))/(2 * SNR))^2



%% Add the two Fts, scaling noise to SNR
FF=noise*FA + FB;
% FF=sqrt(noise * 10 * A^2)*FA + FB; % the factor of 10 because the noise series is 10xlonger than signal pulse
FF = FF * sqrt(sum(sum(abs(FB).^2))/sum(sum(abs(FF).^2))); % scale whole to original power of idat

%ifft to get data:
odat=zeros(nsamps,3);
for ic=1:nchans
    odat(:,ic)=ifft(FF(:,ic));%*A(ic); % for now do not scale by amplitude of original
    
    if plotopt==1
    figure(13);
	amp = max(max(abs(idat)));
    subplot(nchans,2,2*ic); plot([0:dt:samplength-dt],odat(:,ic));
    ylim([-amp  amp])
    subplot(nchans,2,2*ic-1); plot([0:dt:samplength-dt],idat(:,ic));
    ylim([-amp  amp])
    ylabel(sprintf('Component %u',ic));
    end
    
end

%% IF POSS, UNSPLIT AND MEASURE SIGNAL TO NOISE using splitlab method
if nargin > 6
    
dtheta = phi-theta; %I think this is right - this is the clockwise angle by which the true radial is rotated to get to the fast
shiftt = ceil(dT/dt);
dphi = -phi;

de=odat(:,1);
dn=odat(:,2);
%taper
  len  = round(length(de)*.05); %taper length is 3% of total seismogram length
  nn   = 1:len;
  nn2  = (length(de)-len+1):length(de);
  x    = linspace(pi, 2*pi, len);
taper  = 0.5 * (cos(x')+1);
taper2 = fliplr(taper);
de(nn) = de(nn).*taper;     de(nn2) = de(nn2).*taper2;
dn(nn) = dn(nn).*taper;     dn(nn2) = dn(nn2).*taper2;
% filter NB THIS MIGHT NEED TO BE CHANGED DEPENDING ON THE NOISE 
[b,a]  = butter(3, 2*[0.02 0.125]/samprate);
xnf=filtfilt(b,a,dn');
xef=filtfilt(b,a,de');
% unsplit
cr=cosd(-dphi);
sr=sind(-dphi);
rrmat=[cr,sr;-sr,cr];
xf=rrmat(1,:)*[xnf;xef]; xf=xf'; 
xs=rrmat(2,:)*[xnf;xef]; xs=xs';
xf=xf(1:end-shiftt);
xs=xs(shiftt+1:end);

cr=cosd(-dtheta);
sr=sind(-dtheta);
rrmat=[cr,sr;-sr,cr];
xq=rrmat(1,:)*[xf';xs']; xq=xq'; 
xt=rrmat(2,:)*[xf';xs']; xt=xt';
SNRest=max(abs(xq)) / (2*std(xt));
dSNR=abs(SNRest-SNR);
%% if no unsplit parameters available
else 
    dSNR=0;
    for ic=1:nchans
        SNRest(ic)=max(abs(odat(:,ic)))/sqrt(mean(odat(:,ic).^2));
        SNRest(ic)=max(abs(odat(:,ic)))/sqrt(mean(odat(idat(:,ic) < 0.1*max(max(idat)),ic).^2));
        indsig = find(abs(idat(:,ic)) >= var(idat(:,ic))); length(indsig);
        indnoi = find(abs(idat(:,ic)) <= var(idat(:,ic))); length(indnoi); 
        %% till here
        SNRest(ic)=(sum(odat(indsig,ic).^2)/length(indsig)) / (sum(odat(indnoi,ic).^2)/length(indnoi));
        SNRest(ic)=max(abs(odat(:,ic)))/(2*std(odat(indnoi,ic)));
%         SNRest(ic)=sqrt(var(odat(indsig,ic)))/sqrt(var(odat(indnoi,ic)));

    end 
end % if nargin includes unsplit parameters

end % while dSNR too big

%% Degrees of freedom per second
% calculate using first zero-crossing of odat acf
dofps = zeros(nchans,1);
for ic = 1:nchans
    try
    spdof = min(abs(zerof(xcorr(odat(:,ic))) - nsamps))/samprate;
    dofps(ic) = 1./spdof;
    end
end
