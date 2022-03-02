%% 1) Define mouse data
clear all

% data structure:
    % 1) FP raw data
    % 2) EEG raw data
    % 3) EEG sleep score
    % 4) EEG onset clock time (batch I, for alignment)
    % 5) FP onset clock time (batch I, for alignment)
    % 6) time correction for EEG scoring
    
example_mouseID =    {'C:\Users\username\data\FP_data_folder' 'C:\Users\username\data\EEG_data.exp' 'C:\Users\username\data\sleep_score.xlsx' [] [] 0};

mouse = example_mouseID;

%% 2a) Define data streams - Batch I

data = TDTbin2mat(mouse{1});

signal_fs = data.streams.Dv1L.fs; 
signal2 = data.streams.D2B2.data; %hSyn-NE
signal3 = data.streams.Dv1L.data; %FLEX-GCAMP6 in LC
control3 = data.streams.Dv2L.data; %FLEX-GCaMP in LC (violet light)

%% 2b) Define data streams - Batch II

data = TDTbin2mat(mouse{1});

signal_fs = data.streams.Dv1j.fs; 
signal2 = data.streams.Dv2N.data; %hSyn-NE
signal3 = data.streams.D1G2.data; %FLEX-GCAMP6 in LC
control3 = data.streams.D2G2.data; %FLEX-GCaMP in LC (violet light)

% removing FP trace prior to first TTL pulse
TTL_FP = data.epocs.Pu1e.onset;
first_TTL = TTL_FP(1)*signal_fs;
if first_TTL<1
    onset_FP = 1;
end

signal2 = signal2(onset_FP:end);
signal3 = signal3(onset_FP:end);
control3 = control3(onset_FP:end);

%% 3) Normalize and plot FP traces

MeanFilterOrder = 1000; % for smoothing
MeanFilter = ones(MeanFilterOrder,1)/MeanFilterOrder;

fs_signal = 1:length(signal2);
sec_signal = fs_signal/signal_fs;

med_2 = median(signal2);

% deltaF/F
delta2 = ((signal2 - med_2)/med_2)*100;

%normalization of LC-GCaMP to 405nm channel
reg = polyfit(control3(1000*signal_fs:end), signal3(1000*signal_fs:end), 1); 
a = reg(1);
b = reg(2);
controlFit = a.*control3 + b;
controlFit =  filtfilt(MeanFilter,1,double(controlFit));
normDat = (signal3 - controlFit)./controlFit;
delta3 = normDat * 100;

% check fit
figure
a = subplot(4,1,1);
    plot(sec_signal(1000:end), control3(1000:end));
    title('raw control');
b = subplot(4,1,2);
    plot(sec_signal(1000:end), signal3(1000:end));
    title('raw signal');
c = subplot(4,1,3);
    plot(sec_signal(1000:end), signal3(1000:end));
    hold on
    plot(sec_signal(1000:end), controlFit(1000:end));
    title('fitted control');
d = subplot(4,1,4);
    plot(sec_signal(1000:end), delta3(1000:end));
    title('normalized signal');
linkaxes([a,b,c,d],'x');

delta2_filt = filtfilt(MeanFilter,1,double(delta2));
delta3_filt = filtfilt(MeanFilter,1,double(delta3));

% downsampling traces for plotting
ds_delta2_filt = downsample(delta2_filt, 100);
ds_delta3_filt = downsample(delta3_filt, 100);

fs_signal = 1:1:length(delta2_filt);
sec_signal = fs_signal/signal_fs;
ds_sec_signal = downsample(sec_signal, 100); % to remove noise?

% the index 1000:end removes the first second of the recoding for nicer plotting
figure
a = subplot(2,1,1);
    plot(ds_sec_signal(1000:end), ds_delta2_filt(1000:end))
    title('Syn-NE2m');
b = subplot(2,1,2);
    plot(ds_sec_signal(1000:end), ds_delta3_filt(1000:end))
    title('LC-GCaMP');
linkaxes([a,b],'x');


%% 4) loading and plotting EEG and EMG raw data

% Add functions to path
addpath(genpath(['Q:\Personal_folders\Mie\EEG data from NH\EEG toolbox']));

% Import EEG raw data to matlab
Info=loadEXP([mouse{2}],'no');

TimeReldebSec=0; %start extract data from the beginning (first bin)
TimeRelEndSec=inf; %inf to include all data (until last bin)

[Data,Time]=ExtractContinuousData([],Info,[],TimeReldebSec, TimeRelEndSec,[],1);

EMG_rawtrace = Data(1,1:end);
EEG_rawtrace = Data(2,1:end);

sampling_freq = Info.Fs;
EEG_time = (1:length(EEG_rawtrace))/sampling_freq;

% Plot of EEG and EMG traces
figure
a = subplot(2,1,1);
    plot(EEG_time, EMG_rawtrace); 
    xlabel('time (s)');
    ylabel('EMG (V)');
b = subplot(2,1,2);
    plot(EEG_time, EEG_rawtrace); 
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([a,b],'x');


%% 5) open EEG scoring
time_correction = mouse{6};
EEG_sleepscore = xlsread(mouse{3});

wake_onset = rmmissing(EEG_sleepscore(:, 2));
wake_duration = rmmissing(EEG_sleepscore(:, 3)); 

sws_onset = rmmissing(EEG_sleepscore(:, 6)); 
duration_sws = rmmissing(EEG_sleepscore(:, 7)); 

REM_onset = rmmissing(EEG_sleepscore(:, 10)); 
REM_duration = rmmissing(EEG_sleepscore(:, 11)); 

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

wake_binary_vector = zeros([1, (Info.HypnoFiles.Duration)+5]);

for i=1:length(wake_onset) 
    t = wake_onset(i)+1; % +1 to put time 0 as index 1
    d = wake_duration(i)-1; % -1 compensates for adding 1
    wake_binary_vector(t:t+d) = 1;
end

sws_binary_vector = zeros([1, (Info.HypnoFiles.Duration)+5]);

for i=1:length(sws_onset)
    t = sws_onset(i)+1; 
    d = duration_sws(i)-1;
    sws_binary_vector(t:t+d) = 1;
end

REM_binary_vector = zeros([1, (Info.HypnoFiles.Duration)+5]);
for i=1:length(REM_onset) 
    t = REM_onset(i)+1;
    d = REM_duration(i)-1;
    REM_binary_vector(t:t+d) = 1;
end

% Time vector for sleep scoring (1 Hz)
sleepscore_time = 0:length(wake_binary_vector)-1; % Should be same length for wake/sws/REM binary vectors

figure;
a = subplot(2,1,1);
    plot_sleep(EEG_time, EMG_rawtrace, sleepscore_time, wake_binary_vector, sws_binary_vector, REM_binary_vector);
    xlabel('time (s)');
    ylabel('EMG (V)');
b = subplot(2,1,2);
    plot_sleep(EEG_time, EEG_rawtrace, sleepscore_time, wake_binary_vector, sws_binary_vector, REM_binary_vector);
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([a,b],'x');

% 2-column vectors with on- and offsets for each state
wake_periods = [wake_onset wake_onset+wake_duration];
sws_periods = [sws_onset sws_onset+duration_sws];
REM_periods = [REM_onset REM_onset+REM_duration];


%% 6) Dividing wake bouts into microarousals (MA) and wake w/o MA

MA_maxdur = 15; % maximum duration of microarrousal
MA_idx = find(wake_duration < MA_maxdur);
MA_onset = wake_onset(MA_idx);
MA_duration = wake_duration(MA_idx);

MA_binary_vector = zeros([1, (Info.HypnoFiles.Duration)+5]);
for i=1:length(MA_onset)
    t = MA_onset(i)+1;
    d = MA_duration(i)-1;
    MA_binary_vector(t:t+d) = 1;
end

% remove micrarrousal from wake vectors
wake_woMA_onset = wake_onset;
wake_woMA_onset(MA_idx) = [];
wake_woMA_duration = wake_duration;
wake_woMA_duration(MA_idx) = [];
wake_woMA_binary_vector = zeros([1, (Info.HypnoFiles.Duration)+5]); 
for i=1:length(wake_woMA_onset) 
    t = wake_woMA_onset(i)+1;
    d = wake_woMA_duration(i)-1;
    wake_woMA_binary_vector(t:t+d) = 1;
end

% 2-column vectors with on- and offsets for each state
MA_periods = [MA_onset MA_onset+MA_duration];
wake_woMA_periods = [wake_woMA_onset wake_woMA_onset+wake_woMA_duration];


%% 7a) Alingment of EEG recording and FP recording (Batch I)

% This is based on the manually noted time of recording onset
EEG_onset_time = datetime(mouse{4}, 'InputFormat', 'dd-MM-yyyy HH:mm:ss.SSS');
FP_onset_time = datetime(mouse{5}, 'InputFormat', 'dd-MM-yyyy HH:mm:ss.SSS');
time_between_onsets = between(EEG_onset_time, FP_onset_time);

% Number of seconds that should be removed from EEG recording
TTL_EEG_onset = seconds(time(time_between_onsets));

%% 7b) Alingment of EEG recording and FP recording (Batch II)

% TTL pulse from FP
TTL_pulse = Data(3,1:end);
onset_EEG = find(diff(TTL_pulse>1*10^-3));
onset_EEG = onset_EEG(1);

TTL_EEG_onset = onset_EEG/sampling_freq;

%% 7c) Aligning EEG/EMG traces with FP recording

% Cutting EEG/EMG traces leading up to first TTL 
EMG_rawtrace_cut = EMG_rawtrace(round(TTL_EEG_onset*sampling_freq):end);
EEG_rawtrace_cut = EEG_rawtrace(round(TTL_EEG_onset*sampling_freq):end);
EEG_time_cut = (1:length(EEG_rawtrace_cut))/sampling_freq;

wake_binary_vector_cut = wake_binary_vector(round(TTL_EEG_onset+1):end);
sws_binary_vector_cut = sws_binary_vector(round(TTL_EEG_onset+1):end);
REM_binary_vector_cut = REM_binary_vector(round(TTL_EEG_onset+1):end);
    MA_binary_vector_cut = MA_binary_vector(round(TTL_EEG_onset+1):end);
    wake_woMA_binary_vector_cut = wake_woMA_binary_vector(round(TTL_EEG_onset+1):end);

% Align onset, offset, and duration vectors based on TTL
[wake_onset_cut, wake_offset_cut] = binary_to_OnOff(wake_binary_vector_cut);
wake_duration_cut = wake_offset_cut - wake_onset_cut;

[sws_onset_cut, sws_offset_cut] = binary_to_OnOff(sws_binary_vector_cut);
sws_duration_cut = sws_offset_cut - sws_onset_cut;

[REM_onset_cut, REM_offset_cut] = binary_to_OnOff(REM_binary_vector_cut);
REM_duration_cut = REM_offset_cut - REM_onset_cut;

    [MA_onset_cut, MA_offset_cut] = binary_to_OnOff(MA_binary_vector_cut);
    MA_duration_cut = MA_offset_cut - MA_onset_cut;

    [wake_woMA_onset_cut, wake_woMA_offset_cut] = binary_to_OnOff(wake_woMA_binary_vector_cut);
    wake_woMA_duration_cut = wake_woMA_offset_cut - wake_woMA_onset_cut;

% Align period arrays according to TTL
wake_periods_cut = [wake_onset_cut wake_offset_cut];
sws_periods_cut = [sws_onset_cut sws_offset_cut];
REM_periods_cut = [REM_onset_cut REM_offset_cut];
    MA_periods_cut = [MA_onset_cut MA_offset_cut];
    wake_woMA_periods_cut = [wake_woMA_onset_cut wake_woMA_offset_cut];

   
%% 8) Plotting all traces and scorings together

sleepscore_time_cut = 0:length(wake_binary_vector_cut)-1; % should be same length for wake/sws/REM

fig = figure;
a = subplot (4,1,1);
    plot_sleep(ds_sec_signal(1000:end), ds_delta2_filt(1000:end), sleepscore_time_cut, wake_woMA_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut, MA_binary_vector_cut);
    title('Syn-NE2m');
b = subplot(4,1,2);
    plot_sleep(ds_sec_signal(1000:end), ds_delta3_filt(1000:end), sleepscore_time_cut, wake_woMA_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut, MA_binary_vector_cut);
    title('LC-GCaMP');
c = subplot(4,1,3);
    ds_EEG_time = downsample(EEG_time_cut, 10);
    ds_EMG_rawtrace = downsample(EMG_rawtrace_cut, 10);
    plot_sleep(ds_EEG_time, ds_EMG_rawtrace, sleepscore_time_cut, wake_woMA_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut, MA_binary_vector_cut);
    xlabel('time (s)');
    ylabel('EMG (V)');
d = subplot(4,1,4);
    ds_EEG_rawtrace = downsample(EEG_rawtrace_cut, 10);
    plot_sleep(ds_EEG_time, ds_EEG_rawtrace, sleepscore_time_cut, wake_woMA_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut, MA_binary_vector_cut);
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([a,b,c,d],'x');

%% 9) Time points for epoc analysis

Time_points = []; %different transition events (sec) such as EEG-defined transitions, NE trough, NE descend onset prior to REM etc

%% 10) Epoc extraction of LC, NE and EEG power traces

window = 5;

power_bands = {[1, 4], [4, 8], [8, 15], [15, 30]}; 
total_power_band = [0, 30];
frw = 0:0.2:30;
frq = sampling_freq;

time_before = 60; %sec
time_after = 60; %sec 

transition_FP_NE_collector = [];
transition_FP_LC_collector = [];
 
transition_times = arrayfun(@(x) x(1),Time_points); 
transition_times = transition_times(~isnan(transition_times));
transition_spectrogram_collector = [];

for transition_time_number=1:length(transition_times)
    transition_time = transition_times(transition_time_number);
    if (transition_time+time_after)*signal_fs < length(delta2_filt)
        transition_index = round(transition_time*frq);
        transition_index_FP = round(transition_time*signal_fs);
        transition_before_index = round(transition_index-frq*time_before);
        transition_after_index = round(transition_index+frq*time_after);
        transition_before_index_FP = round(transition_index_FP-signal_fs*time_before);
        transition_after_index_FP = round(transition_index_FP+signal_fs*time_after);
        fp_NE_trace = delta2_filt(transition_before_index_FP:transition_after_index_FP);
        fp_LC_trace = delta3_filt(transition_before_index_FP:transition_after_index_FP);
        eeg_transition_trace = EEG_rawtrace_cut(:, transition_before_index:transition_after_index);
        [transition_spectrogram, F, T] = spectrogram(eeg_transition_trace,round(frq*window),[],frw,frq,'yaxis'); % F = frequenciy vector, T=time vector
        transition_spectrogram_collector = cat(3, transition_spectrogram_collector, transition_spectrogram);
        transition_FP_NE_collector = [transition_FP_NE_collector fp_NE_trace'];
        transition_FP_LC_collector = [transition_FP_LC_collector fp_LC_trace'];
    else continue
    end
end

mean_spectrogram = nanmean(log(abs(transition_spectrogram_collector)), 3);

time_spectrogram_zero = T-time_before; 
time_FP = (1:1:length(fp_NE_trace))/signal_fs -time_before;

figure()
a = subplot(6, 1, 1);
    filtered_mean_spectrogram = imgaussfilt(mean_spectrogram, 4);
    imagesc(time_spectrogram_zero, F, filtered_mean_spectrogram); %plot the log spectrum
    set(gca,'YDir', 'normal'); % flip the Y Axis so lower frequencies are at the bottom
    ylim([0, 30]);
    caxis([-6.5, -4.7])
    title('spectrogram');
    colormap(gca, 'parula');
    hold on
    for band_i = 1:length(power_bands)
        plot([-295, -295], power_bands{band_i}, 'LineWidth', 5)
    end
b = subplot(6, 1, 2);
    band_power_collector = [T];
    norm_time = abs(T-(time_before/3*2)); 
    norm_sampling = find(norm_time == min(norm_time));
    normalization_factor = mean(mean(mean_spectrogram(find(F==total_power_band(1)):find(F==total_power_band(2)), 1:norm_sampling)));
    for band_i = 1:length(power_bands)
        power_band = power_bands{band_i};
        power_trace = mean(mean_spectrogram(find(F==power_band(1)):find(F==power_band(2)), :), 1);
        normalized_power_trace = power_trace/-normalization_factor+2;
        band_power_collector = [band_power_collector; normalized_power_trace];
        plot(time_spectrogram_zero, normalized_power_trace)
        hold on
    end
    title('power traces (blue:delta red:theta yellow: alpha (=spindle) purple: beta ');
c = subplot(6, 1, 3);
    sigma = band_power_collector(4,:);
    plot(time_spectrogram_zero, sigma)
    title('sigma');
d = subplot(6, 1, 4);
    NE_mean = mean(transition_FP_NE_collector,2);
    ds_NE_mean = downsample(NE_mean, 100);
    ds_time = downsample(time_FP, 100);
    plot(time_FP, NE_mean )
    title('norepinephrine level');
e = subplot(6, 1, 5);
    LC_mean = mean(transition_FP_LC_collector,2);
    ds_LC_mean = downsample(LC_mean, 100);
    ds_time = downsample(time_FP, 100);
    plot(time_FP, LC_mean )
linkaxes([a,b,c,d,e],'x');


 