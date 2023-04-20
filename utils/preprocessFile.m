function [downsample_signal,t_downsample,is_preterm,name,gestational_age,delivery_age] = preprocessFile(k,myDir)
% preprocess the k-th signal in "myDir", including filtering and
% downsampling.


% Check that k is valid
signalFiles = dir(fullfile(myDir,'*m.mat'));
if k>length(signalFiles)
    error('k exceeds number of files\n');
end

% Get filename
signal_baseFileName = signalFiles(k).name;
name = signal_baseFileName;
signal_fullFileName = fullfile(myDir, signal_baseFileName);

if isequal(signal_baseFileName(end-9),'_')
    header_fullFileName = fullfile(myDir, [signal_baseFileName(1:end-5),'.hea']);
else
    header_fullFileName = fullfile(myDir, [signal_baseFileName(1:end-4),'.hea']);
end

%get days to labor
[is_preterm,gestational_age,delivery_age] = get_preterm_flag(header_fullFileName);

% read signal
load(signal_fullFileName);
fs = 20; %Hz
t0 = 100;

% Bandpass filtering and removing first part of the signal, which sometimes
% contains artifacts.
signal_raw = val([1 3 5 ],t0*fs:end);
t = t0:1/fs:size(val,2)/fs; %s
ehg_sample_signal = signal_raw;
high_frequency_cutoff = 4;
fc_low = high_frequency_cutoff/(fs/2);
[b_low,a_low] = butter(4,fc_low);
fc_high = 0.05/(fs/2);
[b_high,a_high] = butter(4,fc_high,'high');
signal_LP = filtfilt(b_low,a_low,ehg_sample_signal');
signal_filtered = filtfilt(b_high,a_high,signal_LP);

% Downsampling to improve computation speed without losing information.
f_ds = high_frequency_cutoff*2.4;
k = round(fs/f_ds);
downsample_signal = downsample(signal_filtered,k);
t_downsample = downsample(t,k);


end

