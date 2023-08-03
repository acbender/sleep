function [artifact, Sat, EMG, LowV, artifact_pct] = artifact (time, eeg, Fs, SatAmp, EMGcutoff, LVcutoff)

% artifact detector v1

% created April 2022

% uses sliding 5 second window, this can be adjusted

% requires all of the following inputs:
%   time: time-stamps of eeg data in seconds
%   eeg: single-channel eeg data
%   Fs: sampling frequency
%   SatAmp: threshold (in uV) for high voltage cutoff, suggest: 500
%   EMGcutoff: integer multiplier for EMG threshold, suggest: 3
%   LVcutoff: threshold (in arbitrary std units) for low voltage criteria, suggest: 0.005

% outputs:
%   artifact: 1 marks artifact by any criterion at given timestamp, else 0
%   Sat: 1 marks artifact by SatAmp criteria, else 0
%   EMG: 1 marks artifact by EMG criteria, else 0
%   LowV: 1 marks artifact by LowV criteria, else 0
%   artifact_pct: pct of file that was labeled as artifact


%% check
    if length(time)-length(eeg)~=0
        error('error in artifact detector: length of inputs must be equal')
    end

%% filter between 30-50hz for emg
eeg=double(eeg);
[be, ae] = butter(6, [30 50]./(Fs/2)); % bandpass filter
d=eeg;
demg = filtfilt(be,ae,d); %bandpass -- for emg  
std_demg=std(demg);

%% detect artifacts

% using 5 sec sliding window
artifact = zeros([1 length(eeg)]);
Sat = zeros([1 length(eeg)]);
EMG = zeros([1 length(eeg)]);
LowV = zeros([1 length(eeg)]);
win = 5*Fs/2; 
for n = [win+1:length(eeg)-win-1]
    blockEEG = d(n-win:n+win);
    blockEMG = demg(n-win:n+win);
    satEEG = max(abs(blockEEG));
    stdEMG = std(blockEMG);
    stdEEG = std(blockEEG);
    if satEEG>SatAmp
        Sat(n) = 1;
    end
    if stdEMG>(EMGcutoff*std_demg)
        EMG(n) = 1;
    end
    if stdEEG<LVcutoff
        LowV(n) = 1;
    end
    artifact(n) = Sat(n) | EMG(n) | LowV(n);
end
artifact_pct=sum(artifact)/length(artifact);

end
% END

