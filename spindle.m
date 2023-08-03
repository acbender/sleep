function [SPIstart, SPIend, SPImax, SPIpeakfreq, SPImeanfreq, SPIduration, artifact] = spindle (time, eeg, sleep, Fs)

% automatic spindle detector v1
% (excludes artifacts)

% based on Andrillon et al 2011
% created April 2022

% requires all of the following inputs:
%   time: time-stamps of eeg data
%   eeg: single-channel eeg data
%   sleep: sleep-stages (0=wake, 1=N1, 2=N2, 3=N3, 4=REM)
%   Fs: sampling frequency

% outputs:
%   SPIstart: 1 marks start of spindle, else 0
%   SPIend: 1 marks end of spindle, else 0
%   SPImax: 1 marks maximum of spindle amplitude, else 0
%   SPIpeakfreq: frequency at maximum of spindle, else NaN
%   SPImeanfreq: mean frequency of spindle across total duration, else NaN
%   SPIduration: duration of spindle in seconds, else NaN

%   calls function 'artifact.m'


%% check
    if length(time)-length(eeg)~=0 || length(time)-length(sleep)~=0
        error('error in spindle detector: length of inputs must be equal')
    end

    t=time;
    SS2=sleep;
    Fs2=Fs;

%% get artifacts
    SatAmp=400;
    EMGcutoff=3;
    LVcutoff=0.01;
    [artifact, ~, ~, ~, ~] = artifact (time, eeg, Fs, SatAmp, EMGcutoff, LVcutoff);

%% wavelet transform:
    [wt,f]=cwt(eeg,Fs2,'amor','FrequencyLimits',[11 16]); %using Morlet wavelet
    pow=abs(wt); 
    pow2=pow.^2;
    pow=pow2;
    mpow=mean(pow,1);
    MApow=smoothdata(mpow,'movmean',0.1); %smooth w/ moving avg 100 ms window

%get mean/median based on NREM stages only, *excluding artifacts
    meanNREM=median(MApow(SS2>=1 & SS2<=3 & artifact==0)); %using median
    sdNREM=std(MApow(SS2>=1 & SS2<=3 & artifact==0));
    MApowZ=(MApow-meanNREM)./sdNREM;

%% set parameters:

%set amp threshold:
    threshold = 2;
    maxthreshold = 3;
    aboveT = zeros([1 length(t)]);
    aboveT(MApowZ>threshold) = 1;

%set min duration:
    minSPI = round(0.3*Fs2); % sample length for 0.3 sec
    maxSPI = round(3*Fs2);

%% find spindles:
    SPIstart = zeros([1 length(t)]);
    SPIend = zeros([1 length(t)]);
    SPImax = zeros([1 length(t)]);
    SPIpeakfreq = zeros([1 length(t)]);
    SPIstart_time = zeros([1 length(t)]);
    SPIend_time = zeros([1 length(t)]);
    SPIduration = zeros([1 length(t)]);
    SPImeanfreq = zeros([1 length(t)]);

    for n = 2:length(t)-maxSPI
        if SS2(n)>=1 && SS2(n)<=3 && aboveT(n)==1 && aboveT(n-1)==0
            blockMin = aboveT(n:n+minSPI);
            blockMax = aboveT(n:n+maxSPI);
            if sum(blockMin)==length(blockMin) && sum(blockMax)<length(blockMax)
                Slast=find(blockMax==0,1,'first');
                SPIend_i=n+Slast-1;
                [Smax,SPImax_i]=max(MApowZ(n:SPIend_i));
                SPImax_i=n+SPImax_i-1;
                if Smax>maxthreshold && artifact(n)==0
                    SPIstart(n) = 1;
                    SPIstart_time(n)=t(n);
                    SPIend(SPIend_i)=1;
                    SPIend_time(n)=t(SPIend_i);
                    SPImax(SPImax_i)=1;
                    SPIduration(n)=SPIend_time(n)-SPIstart_time(n);
                    SPIpeakfreq(n)=f(find(pow(:,SPImax_i)==max(pow(:,SPImax_i)),1,'first'));
                    [~,maxpow_i]=max(pow(:,n:SPIend_i),[],1);
                    tempf=ones(1,length(maxpow_i));
                    newf=f.*tempf;
                    SPImeanfreq(n)=mean(newf(maxpow_i));
                else
                    SPIstart(n) = 0;
                end
            else
                SPIstart(n) = 0;
            end
        end
    end
    
    %check
    if sum(SPIstart)~=sum(SPIend) || sum(SPIstart)~=sum(SPImax)
        error('Error in spindle detector')
    end

%fix outputs:
SPIpeakfreq(SPIpeakfreq==0)=NaN; %frequency at maximum of spindle, else NaN
SPImeanfreq(SPImeanfreq==0)=NaN; % mean frequency of spindle across total duration, else NaN
SPIduration(SPIduration==0)=NaN; %duration of spindle in seconds, else NaN

end

%END
