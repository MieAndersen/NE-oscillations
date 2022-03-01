%% 1a) Define mouse data

% data structure:
    % 1) FP raw data
    % 2) EEG raw data
    % 3) EEG sleep score
    % 4) 465 channel name
    % 5) 405 channel name (or 560 nm channel)
    % 6) TTL pulse field name
    % 7) laser pulse field name
    % 8) EEG sleepscore time correction
    % 9) Synapse dF/F realtime calculation
    % 10) TTL pulse name name
    % 11) laser pulse onset name
    % 12) laser pulse offset name

example_mouseID = {'C:\Users\username\data\FP_data_folder' 'C:\Users\username\data\EEG_data.exp' 'C:\Users\username\data\sleep_score.xlsx'  'channel 465' 'channel 405' 'TTL pulse' 'laser channel' 0 'Synapse dF/F' 'TTL pulse name' 'laser pulse onset' 'laser pulse offset name'};


mouse = example_mouseID;

%% 2) Load FP data

data = TDTbin2mat(mouse{1}, 'STORE', {mouse{4}, mouse{5}, mouse{9}, mouse{10}, mouse{11}, mouse{12}}); % custom function provided by Tucker Davis Technologies

signal_fs = data.streams.(mouse{4}).fs;

signal_465 = data.streams.(mouse{4}).data; 
signal_405 = data.streams.(mouse{5}).data; 
signal_dF = data.streams.(mouse{9}).data;
TTL_FP = data.epocs.(mouse{6}).onset;
TTL_gap = diff(TTL_FP) > 10 + 1;
if isempty(find(TTL_gap == 1, 1))
    TTL_onset = TTL_FP(1);  % when TTL pulse train is only started once
else 
    TTL_onset = TTL_FP(find(TTL_gap==1)+1); % when TTL pulse train is started more than once
end

first_TTL = TTL_onset(1)*signal_fs;
onset_FP = first_TTL;

laser_on = data.epocs.(mouse{7}).onset-TTL_onset;
laser_on = laser_on( laser_on>=0 );
laser_off = data.epocs.(mouse{7}).offset-TTL_onset;
laser_off = laser_off( laser_off>=0 );

signal_465 = signal_465(onset_FP:end);
signal_405 = signal_405(onset_FP:end);
signal_dF = signal_dF(onset_FP:end);


%% 3) Normalize and plot 

MeanFilterOrder = 1000; 
MeanFilter = ones(MeanFilterOrder,1)/MeanFilterOrder;

fs_signal = 1:1:length(signal_465);
sec_signal = fs_signal/signal_fs;

reg = polyfit(signal_405, signal_465, 1);
a = reg(1);
b = reg(2);
controlFit = a.*signal_405 + b;
controlFit =  filtfilt(MeanFilter,1,double(controlFit));
normDat = (signal_465 - controlFit)./controlFit;
delta_465 = normDat * 100;

% check fit
figure
a = subplot(4,1,1);
plot(sec_signal(1000:end), signal_405(1000:end));
title('raw control');
b = subplot(4,1,2);
plot(sec_signal(1000:end), signal_465(1000:end));
title('raw signal');
c = subplot(4,1,3);
plot(sec_signal(1000:end), signal_465(1000:end));
hold on
plot(sec_signal(1000:end), controlFit(1000:end));
title('fitted control');
d = subplot(4,1,4);
plot(sec_signal(1000:end), delta_465(1000:end));
title('normalized signal');
linkaxes([a,b,c,d],'x');


delta465_filt = filtfilt(MeanFilter,1,double(delta_465));

ds_factor_FP = 100; 
ds_delta465_filt = downsample(delta465_filt, ds_factor_FP);
ds_signal_dF = downsample(signal_dF, ds_factor_FP);
ds_sec_signal = downsample(sec_signal, ds_factor_FP); % for plotting

laser_binary_vector = zeros([1, length(sec_signal)]); 
for i=1:length(laser_on)
    on = laser_on(i)*signal_fs; 
    off = laser_off(i)*signal_fs;
    laser_binary_vector(on:off) = 1;
end

ds_laser_binary_vector = downsample(laser_binary_vector, ds_factor_FP);

% Plot of the three traces above each other (the index 1000:end removes the
% first second of the recoding for nicer plotting)
figure
a = subplot(3,1,1);
plot(ds_sec_signal, detrend(ds_delta465_filt))
title('NE2m');
b = subplot(3,1,2);
plot(ds_sec_signal, ds_signal_dF(1:length(ds_sec_signal)))
title('Synapse dF/F');
c = subplot(3,1,3);
plot(sec_signal, laser_binary_vector)
ylim([-1 2])
title('laser');
linkaxes([a,b,c],'x');

%% 4) loading and plotting EEG and EMG raw data

% Import EEG raw data to matlab
Info=loadEXP([mouse{2}],'no'); % custom function provided by Viewpoint Behavior Technology

TimeReldebSec=0; 
TimeRelEndSec=inf;

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

%% 4b) Interpolate gaps in EEG/EMG 

% Converted recordigns have a ~1s gap in EEG/EMG traces that must be
% removed to run the spectogram function. Here they are interpolated

EEG_NaNs = find(isnan(EEG_rawtrace));

if ~isempty(EEG_NaNs)
    EEG_nans = isnan(EEG_rawtrace); 
    EEG_nans_n = 1:numel(EEG_rawtrace); 
    EEG_intpl = EEG_rawtrace; 
    EEG_intpl(EEG_nans) = interp1(EEG_nans_n(~EEG_nans), EEG_rawtrace(~EEG_nans), EEG_nans_n(EEG_nans));
    EEG_rawtrace = EEG_intpl;
    
    EMG_nans = isnan(EMG_rawtrace); 
    EMG_nans_n = 1:numel(EMG_rawtrace); 
    EMG_intpl = EMG_rawtrace; 
    EMG_intpl(EMG_nans) = interp1(EMG_nans_n(~EMG_nans), EMG_rawtrace(~EMG_nans), EMG_nans_n(EMG_nans));
    EMG_rawtrace = EMG_intpl;
end

% if the last sampling point is a nan and therefore not interpolated, remove this one
if isnan(EEG_rawtrace(end))
    EEG_rawtrace = EEG_rawtrace(1:end-1);
    EMG_rawtrace = EMG_rawtrace(1:end-1);
    EEG_time = EEG_time(1:end-1);
end

%% 5)  Alingment of EEG recording and FP recording

% TTL pulse from FP
TTL_pulse = Data(3,1:end);
onset_EEG = find(diff(TTL_pulse>1*10^-3));
onset_EEG_time = onset_EEG/sampling_freq;
onset_EEG_time_diff = diff(onset_EEG_time);

TTL_gap_EEG = onset_EEG_time_diff > 10;
if isempty(find(TTL_gap_EEG==1, 1))
    onset_EEG = onset_EEG(1);                                               % for #418 skip this if statement and run only this line
else 
    onset_EEG = onset_EEG(find(onset_EEG_time_diff>10)+1);
end

TTL_EEG_onset = onset_EEG./sampling_freq;

%Cutting EEG/EMG traces leading up to first TTL 
EMG_rawtrace_cut = EMG_rawtrace(round(TTL_EEG_onset*sampling_freq):end);
EEG_rawtrace_cut = EEG_rawtrace(round(TTL_EEG_onset*sampling_freq):end);
EEG_time_cut = (1:length(EEG_rawtrace_cut))/sampling_freq;
ds_EEG_time = downsample(EEG_time_cut, 10);

figure
a = subplot(2,1,1);
    plot(ds_sec_signal, detrend(ds_delta465_filt))
    title('NE2m');
b = subplot(2,1,2);
    plot(EEG_time_cut, EMG_rawtrace_cut); 
    xlabel('time (s)');
    ylabel('EMG (V)');
linkaxes([a,b],'x');


%% 6) open EEG scoring

time_correction = mouse{8};
EEG_sleepscore = xlsread(mouse{3});

%Awake
wake_onset = rmmissing(EEG_sleepscore(:, 2));
wake_duration = rmmissing(EEG_sleepscore(:, 3)); 

%Slow-wave sleep
sws_onset = rmmissing(EEG_sleepscore(:, 6)); 
duration_sws = rmmissing(EEG_sleepscore(:, 7)); 

%REM
REM_onset = rmmissing(EEG_sleepscore(:, 10)); 
REM_duration = rmmissing(EEG_sleepscore(:, 11)); 

% Most EEG scorings don't start at time 0 - which shifts the timeline of the
% scoring relative to the EEG/EMG traces - this is corrected for below
if min([wake_onset(1), sws_onset(1), REM_onset(1)]) ~= 0
    EEG_scoring_onset = min([wake_onset(1), sws_onset(1), REM_onset(1)]); 
    wake_onset = wake_onset - EEG_scoring_onset;
    sws_onset = sws_onset - EEG_scoring_onset;
    REM_onset = REM_onset - EEG_scoring_onset;
end

% NB! all EEG/EMG traces are not aligned properly with sleep score - adjust
% time corrections accordingly
wake_onset = wake_onset+time_correction;
sws_onset = sws_onset+time_correction;
REM_onset = REM_onset+time_correction;

% often Info.HypnoFiles.Duration is incorrect and too short - use last bout as end time instead
hypno_duration = max([(wake_onset(end)+wake_duration(end)), (sws_onset(end)+duration_sws(end)), (REM_onset(end)+REM_duration(end))]);

% Create binary vectors for sleep stages
wake_binary_vector = zeros([1, hypno_duration+time_correction]); % vector of zeros matching the length of recording in seconds (+1 for last time interval).
for i=1:length(wake_onset) 
    t = wake_onset(i)+1; % +1 to put time 0 as index 1
    d = wake_duration(i)-1; % -1 compensates for adding 1
    wake_binary_vector(t:t+d) = 1;
end

sws_binary_vector = zeros([1, hypno_duration+time_correction]); % vector of zeros matching the length of recording in seconds (+1 for last time interval).
for i=1:length(sws_onset) 
    t = sws_onset(i)+1; 
    d = duration_sws(i)-1;
    sws_binary_vector(t:t+d) = 1;
end

REM_binary_vector = zeros([1, hypno_duration+time_correction]); % vector of zeros matching the length of recording in seconds (+1 for last time interval).
for i=1:length(REM_onset)
    t = REM_onset(i)+1;
    d = REM_duration(i)-1;
    REM_binary_vector(t:t+d) = 1;
end

% Time vector for sleep scoring (1 Hz)
sleepscore_time = 0:length(wake_binary_vector)-1; % Should be same length for wake/sws/REM binary vectors

% check alignment of scoring and adjust time_correction if necessary
fig = figure;
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


%% 6b) Dividing wake bouts into microarousals (MA) and wake w/o MA

MA_maxdur = 15; % maximum duration of microarrousal
MA_idx = find(wake_duration < MA_maxdur);
MA_onset = wake_onset(MA_idx);
MA_duration = wake_duration(MA_idx);
MA_binary_vector = zeros([1, hypno_duration+time_correction]);
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
wake_woMA_binary_vector = zeros([1, hypno_duration+time_correction]);
for i=1:length(wake_woMA_onset)
    t = wake_woMA_onset(i)+1;
    d = wake_woMA_duration(i)-1;
    wake_woMA_binary_vector(t:t+d) = 1;
end

MA_periods = [MA_onset MA_onset+MA_duration];
wake_woMA_periods = [wake_woMA_onset wake_woMA_onset+wake_woMA_duration];


%% 7) Aligning binary sleep score vectors and on/offsets

% Remove time before TTL from EEG score to align with FP trace
wake_binary_vector_cut = wake_binary_vector(round(TTL_EEG_onset+1):end);
sws_binary_vector_cut = sws_binary_vector(round(TTL_EEG_onset+1):end);
REM_binary_vector_cut = REM_binary_vector(round(TTL_EEG_onset+1):end);

% Align onset, offset, and duration vectors based on TTL
[wake_onset_cut, wake_offset_cut] = binary_to_OnOff(wake_binary_vector_cut);
wake_duration_cut = wake_offset_cut - wake_onset_cut;

[sws_onset_cut, sws_offset_cut] = binary_to_OnOff(sws_binary_vector_cut);
sws_duration_cut = sws_offset_cut - sws_onset_cut;

[REM_onset_cut, REM_offset_cut] = binary_to_OnOff(REM_binary_vector_cut);
REM_duration_cut = REM_offset_cut - REM_onset_cut;

% Align period arrays according to TTL
wake_periods_cut = [wake_onset_cut wake_offset_cut];
sws_periods_cut = [sws_onset_cut sws_offset_cut];
REM_periods_cut = [REM_onset_cut REM_offset_cut];


%% 7b)  Alingment of MA vectors

% Remove time before TTL from EEG score to align with FP trace
MA_binary_vector_cut = MA_binary_vector(round(TTL_EEG_onset+1):end);
wake_woMA_binary_vector_cut = wake_woMA_binary_vector(round(TTL_EEG_onset+1):end);

% Align onset, offset, and duration vectors based on TTL
[MA_onset_cut, MA_offset_cut] = binary_to_OnOff(MA_binary_vector_cut);
MA_duration_cut = MA_offset_cut - MA_onset_cut;

[wake_woMA_onset_cut, wake_woMA_offset_cut] = binary_to_OnOff(wake_woMA_binary_vector_cut);
wake_woMA_duration_cut = wake_woMA_offset_cut - wake_woMA_onset_cut;

MA_periods_cut = [MA_onset_cut MA_offset_cut];
wake_woMA_periods_cut = [wake_woMA_onset_cut wake_woMA_offset_cut];


%% 8) Plotting all traces and scorings together
% plot EEG sleep scoring bouts together with FP data
sleepscore_time_cut = 0:length(wake_binary_vector_cut)-1; % should be same length for wake/sws/REM

fig = figure;
a = subplot(4,1,1);
    plot_sleep(ds_sec_signal, ds_delta465_filt, sleepscore_time_cut, wake_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut, MA_binary_vector_cut);
    title('NE2m');
b = subplot (4,1,2);
    plot_sleep(sec_signal, laser_binary_vector, sleepscore_time_cut, wake_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut, MA_binary_vector_cut);
    ylim([-1 2])
    title('laser');
c = subplot(4,1,3);
    ds_EEG_time = downsample(EEG_time_cut, 10);
    ds_EMG_rawtrace = downsample(EMG_rawtrace_cut, 10);
    plot_sleep(ds_EEG_time, ds_EMG_rawtrace, sleepscore_time_cut, wake_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut, MA_binary_vector_cut);
    xlabel('time (s)');
    ylabel('EMG (V)');
d = subplot(4,1,4);
    ds_EEG_rawtrace = downsample(EEG_rawtrace_cut, 10);
    plot_sleep(ds_EEG_time, ds_EEG_rawtrace, sleepscore_time_cut, wake_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut, MA_binary_vector_cut);
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([a,b,c,d],'x');   

    
%% 9) EEG power spectrum analysis

fs_signal = 1:1:length(delta465_filt);
sec_signal = fs_signal/signal_fs;

EEG_window = 5; %sec. 1 for 30 sec

power_bands = {[1, 4], [4, 8], [8, 15], [15, 30]};
total_power_band = [0, 30];
frw = 0:0.2:30;

Data_EEG = EEG_rawtrace_cut;
frq = sampling_freq;
    
    [transition_spectrogram, F, T] = spectrogram(Data_EEG,round(frq*EEG_window),[],frw,frq,'yaxis');
    mean_spectrogram = log(abs(transition_spectrogram));
    time_spectrogram_zero = T; 
    filtered_mean_spectrogram = imgaussfilt(mean_spectrogram, 4);
    
    specto_fs = length(T)/T(end);
    
figure()
a = subplot(4, 1, 1);
    plot( downsample(sec_signal, 100),downsample(delta465_filt, 100))
 
b = subplot(4, 1, 2);
    imagesc(time_spectrogram_zero, F, filtered_mean_spectrogram); %plot the log spectrum
    set(gca,'YDir', 'normal'); % flip the Y Axis so lower frequencies are at the bottom
    ylim([0, 30]);
    caxis([-6.7, -4])
    colormap(gca, 'parula');
    hold on
    for band_i = 1:length(power_bands)
        plot([-295, -295], power_bands{band_i}, 'LineWidth', 5)
    end
    
c = subplot(4, 1, 3);
    band_power_collector = [T];
    norm_time = abs(T-(100/3*2)); 
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
    
d = subplot(4, 1, 4);
    sigma = band_power_collector(4,:); %theta/delra ratio
    plot(time_spectrogram_zero, sigma)
    title ('sigma')
  
linkaxes([a,b,c,d],'x');
    


%% 10) Plot NE, EEG powerspectrum, sigma and delta/sigma

figure
a = subplot(6,1,1);
    plot(ds_sec_signal, ds_delta465_filt)
    xlabel('time (s)');
    ylabel('dF/F (%)');
    title('NE2m');
b = subplot(6,1,2);
    ds_EEG_rawtrace = downsample(EEG_rawtrace_cut, 10);
    plot(ds_EEG_time, ds_EEG_rawtrace); 
    xlabel('time (s)');
    ylabel('EEG (V)');
c = subplot(6,1,3);
    imagesc(time_spectrogram_zero, F, filtered_mean_spectrogram); %plot the log spectrum
    set(gca,'YDir', 'normal'); % flip the Y Axis so lower frequencies are at the bottom
    ylim([0, 30]);
    caxis([-6.7, -4])
    colormap(gca, 'parula');
    hold on
    for band_i = 1:length(power_bands)
        plot([-295, -295], power_bands{band_i}, 'LineWidth', 5)
    end
d = subplot(6,1,4);
    plot(time_spectrogram_zero, sigma)
    xlabel('time (s)');
    ylabel('power');
    title('sigma');
e = subplot(6,1,5);
    ds_EMG_rawtrace = downsample(EMG_rawtrace_cut, 10);
    plot(ds_EEG_time, ds_EMG_rawtrace); 
    xlabel('time (s)');
    ylabel('EMG (V)');
f = subplot(6,1,6);
    plot(sec_signal, laser_binary_vector)
    ylim([-1 2])
    xlabel('time (s)');
    ylabel('laser power');
    title('laser');
linkaxes([a,b,c,d,e,f],'x');


%% 11) Times of threshold adjusting
% Manually noted times ('h:mm:ss') for chaning threshold adjustments in FP recording.
% time for setting threshold:
    %1: -15
    %2; -10
    %3: -5
    %4: 0
    %5: +5
    %6: No stim

example_thrshldtimes = {'h:mm:ss' 'h:mm:ss' 'h:mm:ss' 'h:mm:ss' 'h:mm:ss' 'h:mm:ss'};

thrshld_start_time = example_thrshldtimes;


%% 12) Defining time intervals for stimulations based accurate time of threshold adjustments

FPrec_UTCstart = data.info.utcStartTime; % UTC time zone 0
FPrec_date = datetime(data.info.date, 'TimeZone', 'Europe/Copenhagen'); % to determine time correction for summer/winter time

if isdst(FPrec_date) == 1
    UTC_correction = 2; %UTC+2
else
    UTC_correction = 1; %UTC+1
end

FPrec_start = duration(FPrec_UTCstart)+hours(UTC_correction); % Start time of FP recording

thrshld_rel_start_time = cellfun(@duration,thrshld_start_time); % start time of threshold relative to recording start
thrshld_abs_start_time = thrshld_rel_start_time+FPrec_start; % absolute start time of threshold relative to recording

for i=1:length(thrshld_rel_start_time)
    if isnan(thrshld_rel_start_time(i))
        thrshld_rel_start_time(i) = duration('11:00:00')+hours(i) - duration(FPrec_start);
    end
end

% defining start and end times of laser stimulation periods
prebase_1_start = seconds(thrshld_rel_start_time(1) - hours(2)) - TTL_onset(1); % 2 h before first stim period
prebase_2_start = seconds(thrshld_rel_start_time(1) - hours(1)) - TTL_onset(1); % 1 h before first stim period
thrs_15_start = seconds(thrshld_rel_start_time(1)) - TTL_onset(1);
thrs_10_start = seconds(thrshld_rel_start_time(2)) - TTL_onset(1);
thrs_5_start = seconds(thrshld_rel_start_time(3)) - TTL_onset(1);
thrs_0_start = seconds(thrshld_rel_start_time(4)) - TTL_onset(1);
thrs_p5_start = seconds(thrshld_rel_start_time(5)) - TTL_onset(1);
postbase_1_start = seconds(thrshld_rel_start_time(6)) - TTL_onset(1); % 1 h after last stim period
postbase_2_start = seconds(thrshld_rel_start_time(6)+ hours(1)) - TTL_onset(1); % 2 h after last stim period

postbase_2_end = seconds(thrshld_rel_start_time(6)+ hours(2)) - TTL_onset(1);

start_time = [prebase_1_start prebase_2_start thrs_15_start thrs_10_start thrs_5_start thrs_0_start thrs_p5_start  postbase_1_start  postbase_2_start];
stop_time = [prebase_2_start thrs_15_start thrs_10_start thrs_5_start thrs_0_start thrs_p5_start  postbase_1_start  postbase_2_start postbase_2_end];

laser_periods = [start_time' stop_time'];

%% Sleep characterization during baseline and stim periods

epoc_time = laser_periods(:,2) - laser_periods(:,1);
collect_PXX = [];
MA_collect_rel_num = [];
MA_collect_abs_num = [];
sws_perc = [];
sws_number = [];
wake_number = [];
sws_boutduration = [];
sigma_mean = [];
rem_perc = [];
power_collect = [];
sigma_data = [];
NE_ampl = [];
perc90sigma = [];
wake_perc =[];

% Find sleep/wake bouts within each baseline/stim period (j)
for j = 1:length(laser_periods)
    t1 = laser_periods(j, 1);
    t2 = laser_periods(j, 2);

    % NREM bouts within periods
    sws_select_onset = sws_onset_cut(find(sws_onset_cut>t1 & sws_onset_cut<t2)); 
    sws_select_offset = sws_offset_cut(find(t1<sws_offset_cut & sws_offset_cut<t2));
    if isempty(sws_select_onset) % in case of no onsets
        sws_select_offset =[];
    elseif  sws_select_onset(1)>sws_select_offset(1)    % if period started during NREM bout, skip the first offset
        sws_select_offset = sws_select_offset(2:end);
    end

    if isempty(sws_select_onset)
      sws_select_onset = [];  
    elseif sws_select_onset(end)>sws_select_offset(end) % if period ended during NREM bout, skip the last onset
        sws_select_onset = sws_select_onset(1:end-1);
    end

    % REM bouts within periods
    rem_select_onset = REM_onset_cut(find(REM_onset_cut>t1 & REM_onset_cut<t2));
    rem_select_offset = REM_offset_cut(find(t1<REM_offset_cut & REM_offset_cut<t2));
    if isempty(rem_select_onset)
        rem_select_offset = [];
    elseif rem_select_onset(1)>rem_select_offset(1)
        rem_select_offset = rem_select_offset(2:end);
    end

    if isempty( rem_select_onset)
       rem_select_onset = [];
    elseif rem_select_onset(end)>rem_select_offset(end)
        rem_select_onset = rem_select_onset(1:end-1);
    end

    % wake excluding MAs
    wake_select_onset = wake_woMA_onset_cut(find(wake_woMA_onset_cut>t1 & wake_woMA_onset_cut<t2));
    wake_select_offset = wake_woMA_offset_cut(find(t1<wake_woMA_offset_cut & wake_woMA_offset_cut<t2));
    if isempty(wake_select_onset)
        wake_select_offset =[];
    elseif isempty(wake_select_offset)
        wake_select_offset = t2;
    elseif wake_select_onset(1)>wake_select_offset(1)
        wake_select_offset = wake_select_offset(2:end);
    end

    if isempty(wake_select_onset)
       wake_select_onset = [];
    elseif   wake_select_onset(end)>wake_select_offset(end)
        wake_select_onset = wake_select_onset(1:end-1);
    end

    % MA bouts within periods
    MA_select_onset = MA_onset_cut(find(MA_onset_cut>t1 & MA_onset_cut<t2));
    MA_select_offset = MA_offset_cut(find(t1<MA_offset_cut & MA_offset_cut<t2));
    if isempty(MA_select_onset)
        MA_select_offset =[];
    elseif MA_select_onset(1)>MA_select_offset(1)
        MA_select_offset = MA_select_offset(2:end);
    end

    if isempty(MA_select_onset)
       MA_select_onset = [];
    elseif   MA_select_onset(end)>MA_select_offset(end)
        MA_select_onset = MA_select_onset(1:end-1);
    end

% Sleep characterization during baseline and stim periods
    sws_diff = sum(sws_select_offset-sws_select_onset);
    rem_diff = sum(rem_select_offset-rem_select_onset);
    wake_diff = sum(wake_select_offset-wake_select_onset); 

    sws_bout_duration = mean(sws_select_offset-sws_select_onset); 
    sws_bout_number = length(sws_select_onset);
    wake_bout_number = length(wake_select_onset);
    sws_perc = [sws_perc sws_diff/epoc_time(j)*100 ];
    wake_perc = [wake_perc wake_diff/epoc_time(j)*100 ];
    rem_perc = [rem_perc rem_diff/epoc_time(j)*100 ];
    ma_select_onset = length((find(MA_onset_cut>t1 & MA_onset_cut<t2)));
    MA_collect_rel_num = [MA_collect_rel_num ma_select_onset/(sws_diff/60)];
    MA_collect_abs_num = [MA_collect_abs_num ma_select_onset];
    wake_number = [wake_number wake_bout_number];
    sws_number = [sws_number sws_bout_number];
    sws_boutduration = [sws_boutduration sws_bout_duration];

% EEG and NE extraction for sigma and amplitude analysis     
    sws_select_onset_idx_EEG = floor(sws_select_onset*sampling_freq);
    sws_select_offset_idx_EEG = floor(sws_select_offset*sampling_freq);
    sws_select_onset_idx_FP = floor(sws_select_onset*signal_fs);
    sws_select_offset_idx_FP = floor(sws_select_offset*signal_fs);
    NREM_data = cell(1, numel(sws_select_onset_idx_EEG));
    NREM_data_collect = [];
    PXX = [];
    
    for i=1:numel(sws_select_onset_idx_EEG)
        NREM_data{i} = EEG_rawtrace_cut(sws_select_onset_idx_EEG(i):sws_select_offset_idx_EEG(i));
        NREM_data_cut = EEG_rawtrace_cut(sws_select_onset_idx_EEG(i):sws_select_offset_idx_EEG(i));
        NREM_data_collect = [NREM_data_collect NREM_data_cut];
        [pxx, f] = pwelch(NREM_data{i}, [], [],[0:0.2:100], sampling_freq); 
        logpxx = 10*log10(pxx);
        FX{i} = f; % frequency vectors 
        PXX(:,i) = logpxx; % PSD estimates for all NREM bouts
    end

    mean_PXX = mean(PXX,2);

    sws_select_onset_idx_sigma = ceil(sws_select_onset/(EEG_window/2));
    sws_select_offset_idx_sigma = floor(sws_select_offset/(EEG_window/2));

    delta = band_power_collector(2,:);
    theta = band_power_collector(3,:);
    beta = band_power_collector(5,:);

    for h=1:numel(sws_select_onset_idx_sigma)
        sigma_data(j,h) = mean(sigma(sws_select_onset_idx_sigma(h):sws_select_offset_idx_sigma(h)));
        perc90sigma(j,h) = prctile(sigma(sws_select_onset_idx_sigma(h):sws_select_offset_idx_sigma(h)),90);
        NE_select = delta465_filt(sws_select_onset_idx_FP(h):sws_select_offset_idx_FP(h));
        perc90 = prctile(NE_select,90);
        perc10 = prctile(NE_select,10);
        NE_ampl(j,h) = perc90-perc10;
    end
    collect_PXX = [collect_PXX mean_PXX];
end

band_estimate = [];
NE_Ampliude_mean = [];
for t = 1:length(laser_periods)
    band_estimate(t,2) = mean( nonzeros(sigma_data(t,:)));
    perc_estimate(t,2) = mean( nonzeros(perc90sigma(t,:)));
    NE_Ampliude_mean(t) = mean( nonzeros(NE_ampl(t,:)));
end

%% Plot NE trace against sleepscore and laser stim

sleepscore_time_cut = 0:length(wake_binary_vector_cut)-1; % should be same length for wake/sws/REM

fig = figure;
a = subplot(2,1,1);
    plot_sleep(ds_sec_signal, ds_delta465_filt, sleepscore_time_cut, wake_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut);
    title('NE2m');
b = subplot (2,1,2);
    plot_sleep(sec_signal, laser_binary_vector, sleepscore_time_cut, wake_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut);
    ylim([-1 2])
    title('laser');
linkaxes([a,b],'x');


%% 8) Find stimulations that occur during NREM

stim_type_on = laser_on;
stim_type_off = laser_off;
stim_period = sws_periods_cut;

stim_during_each = [];
stim_during_any = [];

for i=1:length(stim_type_on) % i is stim idx
    for j=1:length(stim_period) % j is bout idx
    stim_during_each(i,j) = stim_type_on(i) > stim_period(j,1) & stim_type_on(i) < stim_period(j,2);
    stim_during_any(i) = max(stim_during_each(i,1:end));
    end
end

stim_during_any = stim_during_any';

stim_during_period_onset = stim_type_on(logical(stim_during_any));
stim_during_period_offset = stim_type_off(logical(stim_during_any));


%% NE and sigma trace epoc extraction (laser ON)

%Defining laser on/offsets as first and last pulse of pulse train
laser_onset_idx = find( diff(stim_during_period_onset')>0.1)+1;
laser_start = [stim_during_period_onset(1) stim_during_period_onset(laser_onset_idx)'];
laser_offset_idx = find( diff(stim_during_period_offset')>0.1);
laser_stop = [1 stim_during_period_offset(laser_offset_idx)' stim_during_period_offset(end)];

% mean NE at time of laser ON for each threshold
mean_laser_on_465epocs = [];
mean_laser_on_sigmaepocs = [];
for j = 3:7 % only periods with laser stimulations
    t1 = laser_periods(j, 1);
    t2 = laser_periods(j, 2);
    before = 60;
    after = 60;
    laser_on_select = laser_start(find(laser_start>t1 & laser_start<t2));
    laser_on_465epocs = [];
    laser_on_sigmaepocs = [];
    for i = 1:length(laser_on_select)
         time = laser_on_select(i);  
         signal465_epoc = delta465_filt((time - before)*signal_fs:(time + after)*signal_fs);
         sigma_epoc = sigma((time - before)/(EEG_window/2):(time + after)/(EEG_window/2));
         laser_on_465epocs = [laser_on_465epocs signal465_epoc'];
         laser_on_sigmaepocs = [laser_on_sigmaepocs sigma_epoc'];
    end
    mean_laser = mean(laser_on_465epocs,2);
    mean_sigma = mean(laser_on_sigmaepocs,2);
    mean_laser_on_465epocs = [mean_laser_on_465epocs mean_laser];
    mean_laser_on_sigmaepocs = [mean_laser_on_sigmaepocs mean_sigma];
end

fs_signal1 = 1:1:length(mean_laser_on_465epocs);
sec_signal1 = (fs_signal1/signal_fs)-before;
ds_mean_laser_on_465epocs = downsample(mean_laser_on_465epocs, 100);
ds_sec_signal1 = downsample(sec_signal1, 100);
fs_signal2 = 1:1:length(mean_laser_on_sigmaepocs);
sec_signal2 = (fs_signal2*(EEG_window/2))-before;

ds_mean_laser_on_465epocs_norm = ds_mean_laser_on_465epocs - mean(ds_mean_laser_on_465epocs(1:10*signal_fs/100,:),1); % normalized to mean of 10 first seconds of traces (by subtraction)

figure
plot(ds_sec_signal1,ds_mean_laser_on_465epocs(:,1));
hold on
plot(ds_sec_signal1,ds_mean_laser_on_465epocs(:,2))
hold on
plot(ds_sec_signal1,ds_mean_laser_on_465epocs(:,3))
hold on
plot(ds_sec_signal1,ds_mean_laser_on_465epocs(:,4))
hold on
plot(ds_sec_signal1,ds_mean_laser_on_465epocs(:,5))


figure
plot(sec_signal2,mean_laser_on_sigmaepocs(:,1));
hold on
plot(sec_signal2,mean_laser_on_sigmaepocs(:,2))
hold on
plot(sec_signal2,mean_laser_on_sigmaepocs(:,3))
hold on
plot(sec_signal2,mean_laser_on_sigmaepocs(:,4))
hold on
plot(sec_signal2,mean_laser_on_sigmaepocs(:,5))
