%% 1) Define mouse data
clear all

% data structure:
    % 1) EEG raw data
    % 2) EEG sleep score
    % 3) time correction for EEG scoring
    
example_mouseID =    {'C:\Users\username\data\EEG_data.exp' 'C:\Users\username\data\sleep_score.xlsx' 0};

mouse = example_mouseID;

%% 2) loading and plotting EEG and EMG raw data

% Import EEG raw data to matlab
Info=loadEXP([mouse{2}],'no');

TimeReldebSec=0; %start extract data from the beginning (first bin)
TimeRelEndSec=inf; %inf to include all data (until last bin)

[Data,Time]=ExtractContinuousData([],Info,[],TimeReldebSec, TimeRelEndSec,[],1);

EMG_rawtrace = Data(1,1:end);
EEG_rawtrace = Data(2,1:end);
LFP_rawtrace = Data(3,1:end);

%time vector using sampling frequency
sampling_freq = Info.Fs;
EEG_time = (1:length(EEG_rawtrace))/sampling_freq;

% Plot of EEG and EMG traces
figure
a = subplot(3,1,1);
    plot(EEG_time, EMG_rawtrace); 
    xlabel('time (s)');
    ylabel('EMG (V)');
b = subplot(3,1,2);
    plot(EEG_time, EEG_rawtrace); 
    xlabel('time (s)');
    ylabel('EEG (V)');
c = subplot(3,1,3);
    plot(EEG_time, LFP_rawtrace); 
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([a,b,c],'x');


%% 3) open EEG scoring

time_correction = mouse{3}; % NB! If there is a systematic time lag between EEG/EMG traces and scoring adjust for it by seconds here
EEG_sleepscore = xlsread(mouse{2});

wake_onset = rmmissing(EEG_sleepscore(:, 2));
wake_duration = rmmissing(EEG_sleepscore(:, 3));

sws_onset = rmmissing(EEG_sleepscore(:, 6));
duration_sws = rmmissing(EEG_sleepscore(:, 7));

if size(EEG_sleepscore,2) < 9 % in case of no REM bouts
    REM_onset = NaN;
    REM_duration = NaN;
else
    REM_onset = rmmissing(EEG_sleepscore(:, 10));
    REM_duration = rmmissing(EEG_sleepscore(:, 11));
end

% Most EEG scorings don't start at time 0 which shifts the timeline of the
% scoring relative to the EEG/EMG traces - this is corrected for below
if min([wake_onset(1), sws_onset(1), REM_onset(1)]) ~= 0
    EEG_scoring_onset = min([wake_onset(1), sws_onset(1), REM_onset(1)]);
    wake_onset = wake_onset - EEG_scoring_onset;
    sws_onset = sws_onset - EEG_scoring_onset;
    REM_onset = REM_onset - EEG_scoring_onset;
end

% NB! not all EEG/EMG traces are aligned properly with sleep score - adjust
% time corrections accordingly
wake_onset = wake_onset+time_correction; 
sws_onset = sws_onset+time_correction; 
REM_onset = REM_onset+time_correction; 


wake_binary_vector = zeros([1, (Info.HypnoFiles.Duration)+6]);
for i=1:length(wake_onset)
    t = wake_onset(i)+1; % +1 to put time 0 as index 1
    d = wake_duration(i)-1; % -1 compensates for adding 1
    wake_binary_vector(t:t+d) = 1;
end

sws_binary_vector = zeros([1, (Info.HypnoFiles.Duration)+6]);
for i=1:length(sws_onset)
    t = sws_onset(i)+1; 
    d = duration_sws(i)-1;
    sws_binary_vector(t:t+d) = 1;
end

REM_binary_vector = zeros([1, (Info.HypnoFiles.Duration)+6]);
if ~isnan(REM_onset)
    for i=1:length(REM_onset)
        t = REM_onset(i)+1;
        d = REM_duration(i)-1;
        REM_binary_vector(t:t+d) = 1;
    end
end

% Time vector for sleep scoring (1 Hz)
sleepscore_time = 0:length(wake_binary_vector)-1; % Should be same length for wake/sws/REM binary vectors

figure;
h(1) = subplot(2,1,1);
    plot_sleep(EEG_time, EMG_rawtrace, sleepscore_time, wake_binary_vector, sws_binary_vector, REM_binary_vector);
    xlabel('time (s)');
    ylabel('EMG (V)');
h(2) = subplot(2,1,2);
    plot_sleep(EEG_time, EEG_rawtrace, sleepscore_time, wake_binary_vector, sws_binary_vector, REM_binary_vector);
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([h(1),h(2)],'x');

%%2-column vectors with on- and offsets for each state
wake_periods = [wake_onset wake_onset+wake_duration];
sws_periods = [sws_onset sws_onset+duration_sws];
REM_periods = [REM_onset REM_onset+REM_duration];

Data_EEG = EEG_rawtrace;

%% 4) select EEG ad LFP during NREM sleep

sampling_freq = Info.Fs;

sws_EEG_all = [];
sws_LFP_all = [];

%select all incidents of NREM sleep and pool it together.
for i = 1:length(sws_periods)
   period1 = sws_periods(i,:);
   sws_EEG = EEG_rawtrace(period1(1)*sampling_freq:period1(2)*sampling_freq);
   sws_LFP = LFP_rawtrace(period1(1)*sampling_freq:period1(2)*sampling_freq);
   sws_EEG_all = [sws_EEG_all sws_EEG];
   sws_LFP_all = [sws_LFP_all sws_LFP];
end

%% 5) EEG and LFP power spectrograms

SamplingRate =   sampling_freq; %Hz

SampleIdx = 1:length(EEG_rawtrace);
TimeIdx = SampleIdx./SamplingRate;
window = 5; %sec. 1 for 30 sec
slidingwindow = 2.5;
power_bands = {[1, 4], [4, 8], [8, 15], [15, 30]}; 
total_power_band = [0, 30];
frw = 0:0.2:30;

frq = sampling_freq;
    
[transition_spectrogram, F, T] = spectrogram(EEG_rawtrace,round(frq*window),round(slidingwindow*frq),frw,frq,'yaxis');
mean_spectrogram = log(abs(transition_spectrogram));
time_spectrogram_zero = T; 
filtered_mean_spectrogram = imgaussfilt(mean_spectrogram, 4);
specto_fs = length(T)/T(end);

[transition_spectrogram_LFP, F, T1] = spectrogram(LFP_rawtrace,round(frq*window),[],frw,frq,'yaxis');
mean_spectrogram_LFP = log(abs(transition_spectrogram_LFP));
time_spectrogram_zero1 = T1; 
filtered_mean_spectrogram_LFP = imgaussfilt(mean_spectrogram_LFP, 4);
specto_fs1 = length(T1)/T1(end);

    
%% 6) Spindle detection

% filter design
PassBand = [8, 15]; %Hz
[FilterB,FilterA] = butter(3,PassBand*2/SamplingRate);
Filtered_EEG = filtfilt(FilterB,FilterA,LFP_rawtrace);

figure('Name','Filtered EEG');
plot(TimeIdx,LFP_rawtrace);
hold on;
plot(TimeIdx,Filtered_EEG);

figure
plot(TimeIdx,Filtered_EEG);

% rectified (squared) and smoothed and threshold 
Smoothing = 0.25;
AvgOrder = round(Smoothing*SamplingRate); % # of samples to average
MeanWindow = ones(round(AvgOrder),1)./AvgOrder;
SqFiltered_EEG= Filtered_EEG.^2;
AvgSmoothedEEGPower = conv(SqFiltered_EEG,MeanWindow,'same');

figure('Name','Filtered EEG^2 + smoothed');
plot(TimeIdx,SqFiltered_EEG);
hold on;
plot(TimeIdx,AvgSmoothedEEGPower);

StdThreshold = 0.8; 

% detection criteria 
EEG_Threshold = std(AvgSmoothedEEGPower) * StdThreshold;
MinInterval = 0.5; %s 0.5 for spindles
MaxInterval = 10; %s 10 for spindles

% automatically calculated
MinIntervalIdx = 0.5*SamplingRate; %samples
MaxIntervalIdx = 10*SamplingRate; % samples 10 for spindles

% EEG event detection
UpTriggerEvent = TriggerPoints(AvgSmoothedEEGPower,EEG_Threshold,10); %a minimal inter-event interval of 10 samples
DownTriggerEvent = TriggerPointsEnd(AvgSmoothedEEGPower,EEG_Threshold,10);
DetectedInterval = PairEvents(UpTriggerEvent,DownTriggerEvent);
DetectedIntervalDuration = diff(DetectedInterval,1,2);
SelectedInterval = DetectedInterval(DetectedIntervalDuration>MinIntervalIdx & DetectedIntervalDuration<MaxIntervalIdx,:);

figure('Name','detected intervals on AvgSmoothedEEGPower');
plot(TimeIdx,AvgSmoothedEEGPower);
hold on;
for ii=1:size(SelectedInterval,1)
    ttt = SelectedInterval(ii,1):SelectedInterval(ii,2);
    plot(TimeIdx(ttt),AvgSmoothedEEGPower(ttt),'r');
end

figure('Name','detected intervals on AvgSmoothedEEGPower');
plot(TimeIdx,LFP_rawtrace);
hold on;
for ii=1:size(SelectedInterval,1)
    ttt = SelectedInterval(ii,1):SelectedInterval(ii,2);
    plot(TimeIdx(ttt),Data_EEG(ttt),'r');
end

EEG_bandpass = bandpass(LFP_rawtrace,PassBand,512);


%% 7) frequency band power

band_power_collector = [];

norm_time = abs(T-(100/3*2)); 
norm_sampling = find(norm_time == min(norm_time));
normalization_factor = mean(mean(mean_spectrogram(find(F==total_power_band(1)):find(F==total_power_band(2)), 1:norm_sampling)));

figure
for band_i = 1:length(power_bands)
    power_band = power_bands{band_i};
    power_trace = mean(mean_spectrogram(find(F==power_band(1)):find(F==power_band(2)), :), 1);
    normalized_power_trace = power_trace/-normalization_factor+2;
    band_power_collector = [band_power_collector; normalized_power_trace];
    plot(time_spectrogram_zero, normalized_power_trace)
    hold on
end

figure  
a = subplot(5, 1, 1);
    imagesc(time_spectrogram_zero1, F, filtered_mean_spectrogram_LFP); 
    set(gca,'YDir', 'normal'); 
    ylim([0, 15]);
    caxis([-5, -4])
    colormap(gca, 'parula');
b = subplot(5, 1, 2);
    sigma = band_power_collector(3,:);
    delta = band_power_collector(1,:);
    beta = band_power_collector(4,:);
    plot(time_spectrogram_zero, sigma) 
c = subplot(5, 1, 3);
    plot(TimeIdx,LFP_rawtrace);
    hold on;
    for ii=1:size(SelectedInterval,1)
        ttt = SelectedInterval(ii,1):SelectedInterval(ii,2);
        plot(TimeIdx(ttt),Data_EEG(ttt),'r');
    end 
d = subplot(5, 1, 4);
    plot(TimeIdx, EEG_bandpass)
e = subplot(5, 1, 5);
    plot(TimeIdx, EMG_rawtrace)
linkaxes([a,b,c,d,e],'x');


%% 8) Calculate r for sigma veruss spindle occurrences

x1= start_time; % s
x2= end_time; % s

spindle_vector_onset = zeros(length(Data_EEG),1);
for ii=1:size(SelectedInterval,1)
    ttt = SelectedInterval(ii,1);
    spindle_vector_onset(ttt) = 1;
end 
spindle_c_onset = spindle_vector_onset(round(x1*sampling_freq):round(x2* sampling_freq));

sigma_bin = sigma(round(x1*specto_fs):round(x2* specto_fs));

n= 10*specto_fs;
n1= 10*sampling_freq;

sigma_10s_bin = arrayfun(@(i) mean(sigma_bin(i:i+n-1)),1:n:length(sigma_bin)-n+1)'; 
spindleno_10s_bin = arrayfun(@(i) sum(spindle_c_onset(i:i+n1-1)),1:n1:length(spindle_c_onset)-n1+1)'; 

figure
plot(spindleno_10s_bin,sigma_10s_bin);

%% 9) correlate intervals of NREM sleep: sigma and NE

sampling_fs = frq;

x1= start_time; % s
x2= end_time; % s
sigma = band_power_collector(3,:);

downsampl = 100;
ds_N = downsample(delta2_filt(x1*signal_fs:x2*signal_fs),downsampl);
sigma_bin = sigma(round(x1*specto_fs):round(x2* specto_fs));

%Increase number of EEG trace to length of fiber photometry
d2 = interp1(1:numel(sigma_bin),sigma_bin,linspace(1,numel(sigma_bin),numel( ds_N)));
TimeLag = 60;
TimeLagInSamples = TimeLag * signal_fs/downsampl;
Data = d2;
Reference =  ds_N;

[cc1,lags] = xcorr(unity(detrend(Data)), unity(detrend(Reference)),round(TimeLagInSamples),'unbiased');

figure('Name','xcorrs');

plot(lags/(signal_fs/downsampl),cc1);
ylabel('corr. coef.');
xlabel('time lag (s)');
hold on;
plot([0,0],[-1,1],'k'); plot(xlim,[0,0],'k');
ylim([-1,1]);xlim([-TimeLag,TimeLag]);
title('sigma power vs. NE');

collect_corr = [];
time_s = lags2/(signal_fs/downsampl);
collect_corr(:,1) = downsample(time_s, 10);
collect_corr(:,2) = downsample(cc1,10);

