function [SOstart, SOend, SOmax, SOmaxamp, SOfreq, artifact] = slowwave (time, eeg, sleep, Fs, cutoff)

% automatic SO detector v2
% (excludes artifacts)

% based on Nir, Tononi 2011
% created Oct 2022

% requires all of the following inputs:
%   time: time-stamps of eeg data
%   eeg: single-channel eeg data
%   sleep: sleep-stages (0=wake, 1=N1, 2=N2, 3=N3, 4=REM)
%   Fs: sampling frequency
%   cutoff: amplitude threshold in percentile (e.g. '25' will exclude
%   any events with signal below the 25th percentile)

% outputs: ******************************************
%   SOstart: 1 marks start of slow neg half wave, else 0
%   SOend: 1 marks end of slow neg half wave, else 0
%   SOmax: 1 marks negative maximum of slow oscillation half wave, else 0
%   SOfreq: frequency at negative maximum of slow oscillation, else NaN
%   SOmaxamp: envelope amplitude (from Hilbert) at negative maximum of half wave, else NaN

%   calls function 'artifact.m'


%% check
    if length(time)-length(eeg)~=0 || length(time)-length(sleep)~=0
        error('error in detector: length of inputs must be equal')
    end

    t=time;
    SS2=sleep;
    Fs2=Fs;
    eeg=double(eeg);

%% get artifacts
    SatAmp=400;
    EMGcutoff=3;
    LVcutoff=0.01;
    [artifact, ~, ~, ~, ~] = artifact (time, eeg, Fs, SatAmp, EMGcutoff, LVcutoff);

%% filter in SO range

    [be, ae] = butter(2, [0.5 4]./(Fs2/2)); % bandpass filter
    d=eeg;
    eegf = filtfilt(be,ae,d); 
    eegf=eegf-mean(eegf);

    %hilbert to get amplitude
    Hilb = hilbert(eegf);
    SOamp = abs(Hilb);

    % set marker for SO envelope amplitude above Pth percentile (based on
    % NREM stages) using input variable: cutoff
    P = prctile(SOamp(SS2>=1 & SS2<=3),cutoff);
    threshold = P;
    SOaboveT = zeros([1 length(SOamp)]);
    SOaboveT(SOamp>threshold) = 1;

%% find zero crossings

    % upslope:
    zeroUp=islocalmin(cumsum(eegf));
    zeroDown=islocalmax(cumsum(eegf));

%% get neg half-waves

    SOstart = zeros([1 length(t)]);
    SOend = zeros([1 length(t)]);
    SOmax = zeros([1 length(t)]);
    SOfreq = zeros([1 length(t)]);
    SOmaxamp = zeros([1 length(t)]);

    for n=1:length(eeg)-1*Fs2-1
        if SS2(n)>=1 && SS2(n)<=3 && zeroDown(n)
            if sum(zeroUp(n+1:n+1*Fs2+1))>=1
                nextzero=find(zeroUp(n+1:n+1*Fs2+1)==1,1,'first');
                nextzero_i=n+nextzero;
                timediff=time(nextzero_i)-time(n);
                artifact_temp=artifact(n:nextzero_i);
                aboveT=SOaboveT(n:nextzero_i);
                if timediff>=0.25 && timediff<=1.0 && sum(artifact_temp)==0 && sum(aboveT)>0
                    SOstart(n)=1;
                    SOend(nextzero_i)=1;
                    [~,SOmax_i]=min(eegf(n:nextzero_i));
                    SOmax_i=n+SOmax_i;
                    SOmax(SOmax_i)=1;
                    SOmaxamp(SOmax_i)=SOamp(SOmax_i);
                    SOfreq(SOmax_i)=1/(timediff*2);
                    
                end
            end
        end
    end


%fix outputs:
SOfreq(SOfreq==0)=NaN; %frequency at negative maximum of slow oscillation (using a full wave), else NaN
SOmaxamp(SOmaxamp==0)=NaN; %amplitude and negative maximum of negative half wave, else NaN



end

%END
