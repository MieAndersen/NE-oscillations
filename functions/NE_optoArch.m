%% 1a) Define mouse data
clear all
close all

% data structure:
    % 1) FP raw data
    % 2) EEG raw data
    % 3) EEG sleep score
    % 4) 465 channel name
    % 5) 405 channel name
    % 6) TTL pulse channel name
    % 7) laser pulse channel name
    % 8) interval for fitting (polyfit)
    % 9) EEG sleepscore time correction

example_mouseID = {'C:\Users\username\data\FP_data_folder' 'C:\Users\username\data\EEG_data.exp' 'C:\Users\username\data\sleep_score.xlsx' '465 channel' '405 channel name' 'TTL pulse channel name' 'laser pulse channel name' (100:10000) 0};
    
mouse = example_mouseID;

%% 2) Load FP data

data = TDTbin2mat(mouse{1}); % custom function provided by Tucker Davis Technologies

signal_fs = data.streams.(mouse{4}).fs;

signal_465 = data.streams.(mouse{4}).data; %hSyn-NE
signal_405 = data.streams.(mouse{5}).data; %jRGECO

% removing FP trace prior to first TTL pulse
TTL_FP = data.epocs.(mouse{6}).onset;
TTL_gap = diff(TTL_FP) > 5 + 1;
if isempty(find(TTL_gap == 1, 1))
    TTL_onset = TTL_FP(1);  % when TTL pulse train is only started once
else 
    TTL_onset = TTL_FP(find(TTL_gap==1)+1); % when TTL pulse train is started more than once
end

first_TTL = TTL_onset(1)*signal_fs;
onset_FP = first_TTL;

laser_on = data.epocs.(mouse{7}).onset-TTL_onset;
laser_off = data.epocs.(mouse{7}).offset-TTL_onset;

signal_465 = signal_465(onset_FP:end);
signal_405 = signal_405(onset_FP:end);


%% 3) Normalize and plot 

MeanFilterOrder = 1000;
MeanFilter = ones(MeanFilterOrder,1)/MeanFilterOrder;

fs_signal = 1:1:length(signal_465);
sec_signal = fs_signal/signal_fs;

reg = polyfit(signal_405(round(mouse{8}*signal_fs)), signal_465(round(mouse{8}*signal_fs)), 1);
a = reg(1);
b = reg(2);
controlFit = a.*signal_405 + b;
controlFit =  filtfilt(MeanFilter,1,double(controlFit));
normDat = (signal_465 - controlFit)./controlFit;
delta_465 = normDat * 100;


figure % check fit
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
ds_delta465_filt = downsample(delta465_filt, ds_factor_FP); % downsampling traces for plotting
ds_sec_signal = downsample(sec_signal, ds_factor_FP); % for plotting

laser_binary_vector = zeros([1, length(signal_465)]); 
for i=1:length(laser_on)
    on = round(laser_on(i)*signal_fs); 
    off = round(laser_off(i)*signal_fs);
    if off == inf                        % in case recording ended while laser was on
        off = length(signal_465);
    end
    laser_binary_vector(on:off) = 1;
end

ds_laser_binary_vector = downsample(laser_binary_vector, ds_factor_FP);

% Plot of the three traces above each other (the index 1000:end removes the first second of the recoding for nicer plotting)
figure
a = subplot(2,1,1);
plot(ds_sec_signal, ds_delta465_filt)
title('NE2m');
b = subplot(2,1,2);
plot(ds_sec_signal, ds_laser_binary_vector)
ylim([-1 2])
title('laser');
linkaxes([a,b],'x');

%% 4) loading and plotting EEG and EMG raw data

% Import EEG raw data to matlab
Info=loadEXP([mouse{2}],'no'); % custom function provided by Viewpoint Behavior Technology

TimeReldebSec=0; %start extract data from the beginning (first bin)
TimeRelEndSec=inf; %inf to include all data (until last bin)

[Data,Time]=ExtractContinuousData([],Info,[],TimeReldebSec, TimeRelEndSec,[],1);

EMG_rawtrace = Data(1,1:end);
EEG_rawtrace = Data(2,1:end);

%time vector using sampling frequency
sampling_freq = Info.Fs;
EEG_time = (1:length(EEG_rawtrace))/sampling_freq;

% Plot EEG and EMG traces
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

%% 5) open EEG scoring

time_correction = mouse{9};
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


% Most EEG scorings don't start at time 0, which shifts the timeline of the
% scoring relative to the EEG/EMG traces - this is corrected for below
if min([wake_onset(1), sws_onset(1), REM_onset(1)]) ~= 0
    EEG_scoring_onset = min([wake_onset(1), sws_onset(1), REM_onset(1)]); % determines the number of seconds to be subtracted
    wake_onset = wake_onset - EEG_scoring_onset;
    sws_onset = sws_onset - EEG_scoring_onset;
    REM_onset = REM_onset - EEG_scoring_onset;
end

% NB! all EEG/EMG traces are not aligned properly with sleep score - adjust
% time corrections accordingly
wake_onset = wake_onset+time_correction;
sws_onset = sws_onset+time_correction;
REM_onset = REM_onset+time_correction;

% Create binary vectors for sleep stages (frequency = 1 Hz)
wake_binary_vector = zeros([1, (sum([Info.HypnoFiles.Duration]))+time_correction+6]); % vector of zeros matching the length of recording in seconds (+1 for last time interval). Sum is used because of instances where exp was backed up at 18 has two durations
for i=1:length(wake_onset)
    t = wake_onset(i)+1; % +1 to put time 0 as index 1
    d = wake_duration(i)-1; % -1 compensates for adding 1
    wake_binary_vector(t:t+d) = 1;
end

sws_binary_vector = zeros([1, (sum([Info.HypnoFiles.Duration]))+time_correction+6]); % vector of zeros matching the length of recording in seconds  (+1 for last time interval)
for i=1:length(sws_onset)
    t = sws_onset(i)+1; 
    d = duration_sws(i)-1;
    sws_binary_vector(t:t+d) = 1;
end

REM_binary_vector = zeros([1, (sum([Info.HypnoFiles.Duration]))+time_correction+6]); % vector of zeros matching the length of recording in seconds (+1 for last time interval)
for i=1:length(REM_onset)
    t = REM_onset(i)+1;
    d = REM_duration(i)-1;
    REM_binary_vector(t:t+d) = 1;
end

% Time vector for sleep scoring (1 Hz)
sleepscore_time = 0:length(wake_binary_vector)-1;

% check alignment of scoring and adjust time correction if necessary
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


%% 5b) Dividing wake bouts into microarousals (MA) and wake w/o MA

MA_maxdur = 15; % maximum duration of microarrousal
MA_idx = find(wake_duration < MA_maxdur);
MA_onset = wake_onset(MA_idx);
MA_duration = wake_duration(MA_idx);
MA_binary_vector = zeros([1, (sum([Info.HypnoFiles.Duration]))+time_correction+6]);
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
wake_woMA_binary_vector = zeros([1, (sum([Info.HypnoFiles.Duration]))+time_correction+6]);
for i=1:length(wake_woMA_onset) 
    t = wake_woMA_onset(i)+1;
    d = wake_woMA_duration(i)-1;
    wake_woMA_binary_vector(t:t+d) = 1;
end

% 2-column vectors with on- and offsets for each state
MA_periods = [MA_onset MA_onset+MA_duration];
wake_woMA_periods = [wake_woMA_onset wake_woMA_onset+wake_woMA_duration];


%% 6)  Alingment of EEG recording and FP recording

TTL_pulse = Data(3,1:end); % TTL pulse from FP
onset_EEG = find(diff(TTL_pulse>1*10^-3));
onset_EEG_time = onset_EEG/sampling_freq;
onset_EEG_time_diff = diff(onset_EEG_time);

TTL_gap_EEG = onset_EEG_time_diff > 6;
if isempty(find(TTL_gap_EEG==1, 1))
    onset_EEG = onset_EEG(1);
else 
    onset_EEG = onset_EEG(find(onset_EEG_time_diff>5)+1);
end

TTL_EEG_onset = onset_EEG/sampling_freq;

%Cutting EEG/EMG traces leading up to first TTL 
% Remove time before TTL from EEG and EMG raw traces
EMG_rawtrace_cut = EMG_rawtrace(round(TTL_EEG_onset*sampling_freq):end);
EEG_rawtrace_cut = EEG_rawtrace(round(TTL_EEG_onset*sampling_freq):end);
EEG_time_cut = (1:length(EEG_rawtrace_cut))/sampling_freq;

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


%% 6b)  Alingment of MA vectors

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

%% MA during laser stim (for frequency)

laser_ONperiods = [laser_on laser_off];

MA_onset_during_laser = IsInInterval(MA_onset_cut, laser_ONperiods);
n_MA = length(MA_onset_cut);
n_MA_during_laser = sum(MA_onset_during_laser);

laser_sws_intersections = range_intersection(laser_ONperiods',sws_periods_cut'); % custom function by Xavier Xavier
laser_sws_intersect_periods = reshape(laser_sws_intersections,2,length(laser_sws_intersections)/2)'; 
sws_during_laser_dur = laser_sws_intersect_periods(:,2)-laser_sws_intersect_periods(:,1); 
sws_during_laser_cumdur = sum(sws_during_laser_dur); 
total_sws_cumdur = sum(sws_duration_cut);

MA_freq_stim = n_MA_during_laser/(sws_during_laser_cumdur/60); % frequency of MAs when laser is on
MA_freq_between_stim = (n_MA-n_MA_during_laser)/((total_sws_cumdur-sws_during_laser_cumdur)/60);  % frequency of MAs when laser is off

%% 7) Plotting all traces and scorings together

% plot EEG sleep scoring bouts together with FP data

% Time vector for sleep scoring (1 Hz)

sleepscore_time_cut = 0:length(wake_binary_vector_cut)-1; % should be same length for wake/sws/REM

fig = figure;
a = subplot(4,1,1);
    plot_sleep(ds_sec_signal, ds_delta465_filt, sleepscore_time_cut, wake_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut);
    title('NE2m');
b = subplot (4,1,2);
    plot_sleep(ds_sec_signal, ds_laser_binary_vector, sleepscore_time_cut, wake_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut);
    ylim([-1 2])
    title('laser');
c = subplot(4,1,3);
    ds_EEG_time = downsample(EEG_time_cut, 10);
    ds_EMG_rawtrace = downsample(EMG_rawtrace_cut, 10);
    plot_sleep(ds_EEG_time, ds_EMG_rawtrace, sleepscore_time_cut, wake_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut);
    xlabel('time (s)');
    ylabel('EMG (V)');
d = subplot(4,1,4);
    ds_EEG_rawtrace = downsample(EEG_rawtrace_cut, 10);
    plot_sleep(ds_EEG_time, ds_EEG_rawtrace, sleepscore_time_cut, wake_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut);
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([a,b,c,d],'x');


%% Overall sleep charactertization
% Total duration of sleep stages from first laser stimulation and onwards
% During which state is laser stimulations initiated during

stim_start = laser_on(1);
recording_end = max([wake_woMA_offset_cut(end), sws_offset_cut(end), REM_offset_cut(end), MA_offset_cut(end)]);

% total time (in sec) from first stim until end of recording
poststim_totaltime = recording_end-laser_on(1); 

% bout duration after first stim
poststim_sws_dur = sws_duration_cut(sws_onset_cut>stim_start);
poststim_wake_woMA_dur = wake_woMA_duration_cut(wake_woMA_onset_cut>stim_start);
poststim_REM_dur = REM_duration_cut(REM_onset_cut>stim_start);
poststim_MA_dur = MA_duration_cut(MA_onset_cut>stim_start);

% Creating one vector with different behaviors represented by unique
% numbers (1=wake, 4=sws, 9=REM, 15=MA) at frequency 1Hz
boutscore_vector = zeros([1, sum(Info.HypnoFiles.Duration)+1]);

    % Here using the aligned "cut" vectors
    for i=1:length(wake_woMA_onset_cut)
        t = wake_woMA_onset_cut(i)+1;
        d = wake_woMA_duration_cut(i)-1;
        boutscore_vector(t:t+d) = 1; % wake=1
    end

    for i=1:length(sws_onset_cut)
        t = sws_onset_cut(i)+1;
        d = sws_duration_cut(i)-1;
        boutscore_vector(t:t+d) = 4; % sws=4
    end

    if ~isnan(REM_onset)
        for i=1:length(REM_onset_cut)
            t = REM_onset_cut(i)+1;
            d = REM_duration_cut(i)-1;
            boutscore_vector(t:t+d) = 9; %REM=9
        end
    end

    for i=1:length(MA_onset_cut)
        t = MA_onset_cut(i)+1;
        d = MA_duration_cut(i)-1;
        boutscore_vector(t:t+d) = 15; %MA=15
    end

% trim and include the first bout (the one during which first stimulation starts)
if boutscore_vector(round(stim_start)) == 1 % i.e. if first stim starts during wake bout
    bouttime_before_stim = stim_start-wake_woMA_onset_cut(find(wake_woMA_onset_cut>stim_start,1)-1);
    poststim_wake_woMA_dur = [wake_woMA_duration_cut(find(wake_woMA_onset_cut>stim_start,1)-1)-bouttime_before_stim; poststim_wake_woMA_dur];
elseif boutscore_vector(round(stim_start)) == 4 % i.e. if first stim starts during sws bout
    bouttime_before_stim = stim_start-sws_onset_cut(find(sws_onset_cut>stim_start,1)-1);
    poststim_sws_dur = [sws_duration_cut(find(sws_onset_cut>stim_start,1)-1)-bouttime_before_stim; poststim_sws_dur];
elseif boutscore_vector(round(stim_start)) == 9 % i.e. if first stim starts during REM bout
    bouttime_before_stim = stim_start-REM_onset_cut(find(REM_onset_cut>stim_start,1)-1);
    poststim_REM_dur = [REM_duration_cut(find(REM_onset_cut>stim_start,1)-1)-bouttime_before_stim; poststim_REM_dur];
end

% total duration (% of total recording from first stim)
poststim_totaldur_sws = sum(poststim_sws_dur)/poststim_totaltime*100;
poststim_totaldur_wake = sum(poststim_wake_woMA_dur)/poststim_totaltime*100;
poststim_totaldur_REM = sum(poststim_REM_dur)/poststim_totaltime*100;
poststim_totaldur_MA = sum(poststim_MA_dur)/poststim_totaltime*100;

% which state is animal in at onset of laser stimulations (1=wake, 4=sws, 9=REM, 15=MA)
laser_on_boutscore = boutscore_vector(round(laser_on));

total_laser_number = length(laser_on_boutscore);
laser_on_wake = sum(laser_on_boutscore==1)/total_laser_number*100;
laser_on_sws = sum(laser_on_boutscore==4)/total_laser_number*100;
laser_on_REM = sum(laser_on_boutscore==9)/total_laser_number*100;
laser_on_MA = sum(laser_on_boutscore==15)/total_laser_number*100;

   
%% 8) Find stimulations that occur during NREM

stim_type_on = laser_on;
stim_type_off = laser_off;
stim_period = sws_periods_cut;          % select brain state stims should occur during
stim_strigger_trans = REM_onset_cut;    % select brain state onset stimulation should result in (for next section)

stim_during_each = [];
stim_during_any = [];

for i=1:length(stim_type_on) % i is stim idx
    for j=1:length(stim_period) % j is bout idx
    stim_during_each(i,j) = stim_type_on(i) > stim_period(j,1) & stim_type_on(i) < stim_period(j,2); 
    stim_during_any(i) = max(stim_during_each(i,1:end));
    end
end

stim_during_any = stim_during_any'; % logical vecotr with 1 indicating stims that occur during selected period

stim_during_period_onset = stim_type_on(logical(stim_during_any));
stim_during_period_offset = stim_type_off(logical(stim_during_any));

pct_totalstim = length(stim_during_period_onset)/length(stim_type_on)*100;

%% 9) Find stimulations during which animal transitions to REM/NREM sleep

max_time2transition = 120; % maximum time (s) from stim

stim_trigger_trans = logical(zeros([1,length(stim_during_period_onset)])');
stim_trigger_NO_trans = logical(zeros([1,length(stim_during_period_onset)])');
times2transitions = [];
for i=1:length(stim_during_period_onset)
    next_transition = stim_strigger_trans(find(stim_strigger_trans > stim_during_period_onset(i), 1));  % this returns the first transition after the stim
    time2transition = next_transition - stim_during_period_onset(i);
    if time2transition < max_time2transition 
        stim_trigger_trans(i) = 1;
    end

    if isempty(time2transition)
        time2transition = NaN;
    end
    times2transitions(i) = time2transition;
end

stim_trigger_trans_onset = stim_during_period_onset(stim_trigger_trans);
stim_trigger_NO_trans = ~stim_trigger_trans;
stim_trigger_NO_wakeup_onset =  stim_during_period_onset(stim_trigger_NO_trans); 

n_REMbouts = length(REM_onset);
n_stim_during_period = length(stim_during_period_onset);
n_stim_trigger_REM = length(stim_trigger_trans_onset);

mean_time_to_REM = mean(times2transitions(stim_trigger_trans));

%% 10) NE trace epoc extraction (laser ON/OFF)

before = 60;
after = 90;
downsample_prism = 70;

on_times = stim_during_period_onset; 
off_times = stim_during_period_offset; 

laser_on_465epocs = [];
for i = 1:length(on_times)
     time = on_times(i);  
     signal465_epoc = delta465_filt((time - before)*signal_fs:(time + after)*signal_fs);
     laser_on_465epocs(:,i) = signal465_epoc;
end

laser_off_465epocs = [];
for i = 1:length(off_times)
     time = off_times(i);  
     signal465_epoc = delta465_filt((time - before)*signal_fs:(time + after)*signal_fs);
     laser_off_465epocs(:,i) = signal465_epoc;
end

fs_signal = 1:1:length(signal465_epoc);
sec_signal = (fs_signal/signal_fs)-before;

laser_on_465epocs_mean_ds = downsample(mean(laser_on_465epocs,2), downsample_prism);
laser_off_465epocs_mean_ds = downsample(mean(laser_off_465epocs,2), downsample_prism);

sec_signal_ds = downsample(sec_signal, downsample_prism)';

figure
a = subplot(2,1,1);
    plot(sec_signal_ds, laser_on_465epocs_mean_ds);    
    title('laser on - mean NE')
b = subplot(2,1,2);
    plot(sec_signal_ds, laser_off_465epocs_mean_ds);    
    title('laser off - mean NE')
linkaxes([a,b],'x');

pre_interval = (-30:0);

prebase_start = 50; 
postbase_start = 10; 
postbase_end = 60; 

%summary data of laser on
laser_on_465epocs_mean = mean(laser_on_465epocs,2);
NE_baseline_pre_on = mean(laser_on_465epocs_mean(round((before-prebase_start)*signal_fs:before*signal_fs)));
NE_baseline_post_on = prctile(laser_on_465epocs_mean(round((before+postbase_start)*signal_fs:(before+postbase_end)*signal_fs)),5);
NE_diff_prepost_on = NE_baseline_post_on-NE_baseline_pre_on;

%summary data of laser off
laser_off_465epocs_mean = mean(laser_off_465epocs,2);
NE_baseline_pre_off = mean(laser_off_465epocs_mean(round((before-prebase_start)*signal_fs:before*signal_fs)));
NE_baseline_post_off = prctile(laser_off_465epocs_mean(round((before+postbase_start)*signal_fs:(before+postbase_end)*signal_fs)),95);
NE_diff_prepost_off = NE_baseline_post_off-NE_baseline_pre_off;

% mean NE during full laser on
after_full = before+120;
mean_full_trace = [mean(laser_on_465epocs,2);mean(laser_off_465epocs,2)];
mean_NE_laser_on = mean(mean_full_trace(round((before+postbase_start)*signal_fs:after_full*signal_fs)));

figure
a = subplot(2,1,1);
    plot(sec_signal_ds, downsample(laser_on_465epocs, downsample_prism));
    hold on
    plot(sec_signal_ds, laser_on_465epocs_mean_ds,'LineWidth',1,'Color',[0 0 0]);
    title('laser on')
b = subplot(2,1,2);
    plot(sec_signal_ds, downsample(laser_off_465epocs, downsample_prism));
    hold on
    plot(sec_signal_ds, laser_off_465epocs_mean_ds,'LineWidth',1,'Color',[0 0 0]);
    title('laser off')
linkaxes([a,b],'x');


%% 11) EEG power spectrum analysis
fs_signal = 1:1:length(delta465_filt);
sec_signal = fs_signal/signal_fs;

window = 5; %sec. 1 for 30 sec

power_bands = {[1, 4], [4, 8], [8, 15], [15, 30]};
total_power_band = [0, 30];
frw = 0:0.2:30;

Data_EEG = EEG_rawtrace_cut;
frq = sampling_freq;
    
[transition_spectrogram, F, T] = spectrogram(Data_EEG,round(frq*window),[],frw,frq,'yaxis');
mean_spectrogram = log(abs(transition_spectrogram));
time_spectrogram_zero = T; 
filtered_mean_spectrogram = imgaussfilt(mean_spectrogram, 4);

specto_fs = length(T)/T(end);

figure()
a = subplot(4, 1, 1);
    plot( downsample(sec_signal, 100),downsample(delta465_filt, 100))
 
b = subplot(4, 1, 2);
    imagesc(time_spectrogram_zero, F, filtered_mean_spectrogram); 
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
    normalization_factor = mean(mean(mean_spectrogram(find(F==total_power_band(1)):find(F==total_power_band(2)), 1:find(T==round(1)))));
    for band_i = 1:length(power_bands)
        power_band = power_bands{band_i};
        power_trace = mean(mean_spectrogram(find(F==power_band(1)):find(F==power_band(2)), :), 1);
        normalized_power_trace = power_trace;
        band_power_collector = [band_power_collector; normalized_power_trace];
        plot(time_spectrogram_zero, normalized_power_trace)
        hold on
    end
    
d = subplot(4, 1, 4);
    sigma = band_power_collector(4,:);
    plot(time_spectrogram_zero, sigma)
    title ('sigma')
  
linkaxes([a,b,c,d],'x');
    

%% 12) Plot NE, EEG powerspectrum, sigma and delta/sigma

figure
a = subplot(5,1,1);
    plot(ds_sec_signal, ds_delta465_filt)
    xlabel('time (s)');
    ylabel('dF/F (%)');
    title('NE2m');
b = subplot(5,1,2);
    plot(ds_sec_signal, ds_laser_binary_vector)
    ylim([-1 2])
    xlabel('time (s)');
    ylabel('laser power');
    title('laser');
c = subplot(5,1,3);
	filtered_mean_spectrogram); %plot the log spectrum
    set(gca,'YDir', 'normal');
    ylim([0, 30]);
    caxis([-6.7, -4])
    colormap(gca, 'parula');
    hold on
    for band_i = 1:length(power_bands)
        plot([-295, -295], power_bands{band_i}, 'LineWidth', 5)
    end
d = subplot(5,1,4);
    plot(time_spectrogram_zero, sigma)
    xlabel('time (s)');
    ylabel('power');
    title('sigma');
e = subplot(5,1,5);
    ds_EEG_rawtrace = downsample(EEG_rawtrace_cut, 10);
    plot_sleep(ds_EEG_time, ds_EEG_rawtrace, sleepscore_time_cut, wake_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut);
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([a,b,c,d,e],'x');

   
   
%% 13a) EEG sigma epoc extraction

Time_points = stim_during_period_onset;

time_before = 60; %sec
time_after = 95; %sec 
 
analysis_window = 5; %sec. 1 for 30 sec

power_bands = {[1, 4], [4, 8], [8, 15], [15, 30] [6 8]}; % delta was 0.2 before
%power_bands = {[0.2, 4], [4, 7], [7, 15], [15, 30] [7 10]};
total_power_band = [0, 30];
frw = 0:0.2:30;
frq = sampling_freq;

ds_NE_mean = [];
sigma_epocs = [];
    
sound_FP_NE_collector = [];

sound_spectrogram_collector = [];

for sound_time_number=1:length(Time_points)
    sound_time = Time_points(sound_time_number);
    if (sound_time+time_after)*signal_fs < length(delta465_filt)
        sound_index = round(sound_time*frq);
        sound_index_FP = round(sound_time*signal_fs);
        sound_before_index = round(sound_index-frq*time_before);
        sound_after_index = round(sound_index+frq*time_after);
        sound_before_index_FP = round(sound_index_FP-signal_fs*time_before);
        sound_after_index_FP = round(sound_index_FP+signal_fs*time_after);
        fp_NE_trace = delta465_filt(sound_before_index_FP:sound_after_index_FP);
        eeg_sound_trace = EEG_rawtrace_cut(:, sound_before_index:sound_after_index);
        [sound_spectrogram, F, T] = spectrogram(eeg_sound_trace,round(frq*analysis_window),[],frw,frq,'yaxis'); % F = frequenciy vector, T=time vector
        sigma_power_trace = mean(sound_spectrogram(find(F==8):find(F==15), :), 1);
        normalized_sigma_power = sigma_power_trace/-normalization_factor+2;
        sigma_epocs = [sigma_epocs  normalized_sigma_power'];
        sound_spectrogram_collector = cat(3, sound_spectrogram_collector, sound_spectrogram);
        sound_FP_NE_collector = [sound_FP_NE_collector fp_NE_trace'];
    else continue
    end
end

mean_spectrogram = nanmean(log(abs(sound_spectrogram_collector)), 3);

time_spectrogram_zero = T-time_before; % to get REM onset at 0 instead of 300
time_FP = (1:1:length(fp_NE_trace))/signal_fs -time_before;

figure()
a = subplot(4, 1, 1);
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
b = subplot(4, 1, 2);
    band_power_collector = [T];
    norm_time = abs(T-(time_before/3*2)); 
    norm_sampling = find(norm_time == min(norm_time));
    normalization_factor = mean(mean(mean_spectrogram(find(F==total_power_band(1)):find(F==total_power_band(2)), 1:norm_sampling)));
    for band_i = 1:length(power_bands)
        power_band = power_bands{band_i};
        sigma_power_trace = mean(mean_spectrogram(find(F==power_band(1)):find(F==power_band(2)), :), 1);
        normalized_sigma_power = sigma_power_trace/-normalization_factor+2;
        band_power_collector = [band_power_collector; normalized_sigma_power]; % matrix of timevector and aligned frequency band power traces
        plot(time_spectrogram_zero, normalized_sigma_power)
        hold on
    end
     title('power traces (blue:delta red:theta yellow: alpha (=spindle) purple: beta ');
c = subplot(4, 1, 3);
    sigma_mean = band_power_collector(4,:); %sigma
    plot(time_spectrogram_zero, sigma_mean)
    title('sigma');
d = subplot(4, 1, 4);
    NE_mean = mean(sound_FP_NE_collector,2);
    ds_NE_mean = downsample(NE_mean, 100);
    ds_time = downsample(time_FP, 100);
    plot(time_FP, NE_mean )
    title('norepinephrine level');
linkaxes([a,b,c,d],'x');
 
sigma_mean_60pre_laser = mean(sigma_mean(round((time_before-55)*specto_fs:time_before*specto_fs))); % mean sigma 55 s pre laser
sigma_max_60post_laser = prctile(sigma_mean(round(time_before*specto_fs:(time_before+60)*specto_fs)),95); % top 95 percentile of sigma trace 60 s post laser

delta_mean = band_power_collector(2,:);
delta_mean_60pre_laser = mean(delta_mean(round((time_before-55)*specto_fs:time_before*specto_fs))); % mean sigma 55 s pre laser
delta_max_60post_laser = prctile(delta_mean(round(time_before*specto_fs:(time_before+60)*specto_fs)),95); % top 95 percentile of sigma trace 60 s post laser

theta_mean = band_power_collector(3,:);
theta_mean_60pre_laser = mean(theta_mean(round((time_before-55)*specto_fs:time_before*specto_fs))); % mean sigma 55 s pre laser
theta_max_60post_laser = prctile(theta_mean(round(time_before*specto_fs:(time_before+60)*specto_fs)),95); % top 95 percentile of sigma trace 60 s post laser

%% 13b) theta during REM

theta_data = [];

REM_onset_idx_theta = ceil(REM_onset_cut/(analysis_window/2)); 
REM_offset_idx_theta = floor(REM_offset_cut/(analysis_window/2));

for h=1:numel(REM_onset_idx_theta)
    theta_data(h) = mean(theta(REM_onset_idx_theta(h):REM_offset_idx_theta(h)));
end

REM_theta_estimate = mean(nonzeros(theta_data));

%% 14) Re-classify MAs as NREM sleep (for NE oscillation rate)

% Creating one vector with different behaviors represented by unique
% numbers (1=wake, 4=sws, 9=REM, 15=MA) at frequency 1Hz
boutscore_vector = zeros([1, sum(Info.HypnoFiles.Duration)+time_correction+6]);

for i=1:length(wake_woMA_onset)
    t = wake_woMA_onset(i)+1;
    d = wake_woMA_duration(i)-1;
    boutscore_vector(t:t+d) = 1; % wake=1
end

for i=1:length(sws_onset)
    t = sws_onset(i)+1;
    d = duration_sws(i)-1;
    boutscore_vector(t:t+d) = 4; % sws=4
end

if ~isnan(REM_onset)
    for i=1:length(REM_onset)
        t = REM_onset(i)+1;
        d = REM_duration(i)-1;
        boutscore_vector(t:t+d) = 9; %REM=9
    end
end

for i=1:length(MA_onset)
    t = MA_onset(i)+1;
    d = MA_duration(i)-1;
    boutscore_vector(t:t+d) = 15; %MA=15
end

% re-classify MA as NREM
NREMinclMA_binary_vector = boutscore_vector==4 | boutscore_vector==15;

% align new NREM vector
NREMinclMA_binary_vector_cut = NREMinclMA_binary_vector(round(TTL_EEG_onset+1):end);
[NREMinclMA_onset_cut, NREMinclMA_offset_cut] = binary_to_OnOff(NREMinclMA_binary_vector_cut);
NREMinclMA_duration = NREMinclMA_offset_cut-NREMinclMA_onset_cut;
NREMinclMA_periods_cut = [NREMinclMA_onset_cut NREMinclMA_offset_cut];



%% 15) NE osciallation rate

% Divide NREM into laser on/off
end_time = min(length(signal_465)/signal_fs, length(EEG_rawtrace)/sampling_freq); %
laser_ONperiods = [laser_on laser_off];
laser_OFFperiods = laser_ONperiods-300; % use 2 min periods of laser off to counterbalance laser on periods

timetrace = (1:length(delta465_filt))/signal_fs;
pk_times = [];

min_pkDist = 10; % use 10
min_pkProm =  1.5; % use 1.5 
min_pkWdth = 0; % use 0
MeanFilterOrder_i = 1000;
MeanFilter_i = ones(MeanFilterOrder_i,1)/MeanFilterOrder_i;


fs_signal = 1:1:length(delta465_filt);
sec_signal = fs_signal/signal_fs;
if sec_signal(end) < NREMinclMA_periods_cut(end) % if last NREM bout goes further than FP recording
    NREMinclMA_periods_cut(end) = sec_signal(end); 
end

for i = 1:length(NREMinclMA_periods_cut)
    period_duration_i = NREMinclMA_periods_cut(i,2)-NREMinclMA_periods_cut(i,1);
    if period_duration_i < 60 % periods shorter than 60 s are excluded from analysis
        continue
    end
    period_i = NREMinclMA_periods_cut(i,:)*signal_fs;
    if period_i(1) == 0 % if first NREM bout start at time=0, use idx=1
        period_i(1) = 1;
    end
    NEtrace_i = delta465_filt(round((period_i(1):period_i(2))));
    timetrace_i = timetrace(round((period_i(1):period_i(2))));
    [pks, pklocs, w, p] = findpeaks(NEtrace_i, timetrace_i, 'MinPeakDistance', min_pkDist, 'MinPeakWidth', min_pkWdth, 'MinPeakProminence',min_pkProm); 
    pk_times = [pk_times pklocs];
    figure
    set(gcf, 'Position',  [100, 300, 1500, 250])
    findpeaks(NEtrace_i, timetrace_i, 'MinPeakDistance', min_pkDist,'MinPeakWidth', min_pkWdth, 'MinPeakProminence',min_pkProm, 'Annotate','extents','WidthReference','halfheight');

end

laserON_pk_times = [];
laserON_pks_n = [];
for i = 1:length(laser_ONperiods)
    laserON_pks_n_i = sum(pk_times > laser_ONperiods(i,1) & pk_times <laser_ONperiods(i,2));
    laserON_pks_n = [laserON_pks_n laserON_pks_n_i];
    laserON_pks_idx_i = find(pk_times > laser_ONperiods(i,1) & pk_times <laser_ONperiods(i,2));
    laserON_pk_times = [laserON_pk_times pk_times(laserON_pks_idx_i)];
end

laserOFF_pk_times = [];
laserOFF_pks_n = [];
for i = 1:length(laser_OFFperiods)
    laserOFF_pks_n_i = sum(pk_times > laser_OFFperiods(i,1) & pk_times <laser_OFFperiods(i,2));
    laserOFF_pks_n = [laserOFF_pks_n laserOFF_pks_n_i];
    laserOFF_pks_idx_i = find(pk_times > laser_OFFperiods(i,1) & pk_times <laser_OFFperiods(i,2));
    laserOFF_pk_times = [laserOFF_pk_times pk_times(laserOFF_pks_idx_i)];
end

% check detection in ON/OFF periods
fig = figure;
a = subplot(2,1,1);
    plot(sec_signal, delta465_filt);
    hold on
    plot (laserON_pk_times, delta465_filt(laserON_pk_times*signal_fs), 'r*')
    plot (laserOFF_pk_times, delta465_filt(round(laserOFF_pk_times*signal_fs)), 'b*')
    title('NE2m');
b = subplot (2,1,2);
    plot(sec_signal, laser_binary_vector);
    ylim([-1 2])
    title('laser');
linkaxes([a,b],'x');


laserON_pks_total = sum(laserON_pks_n);
laserOFF_pks_total = sum(laserOFF_pks_n);

% total duration of sws during laser ON periods
laserON_sws_intersections = range_intersection(laser_ONperiods',NREMinclMA_periods_cut'); 
laserON_sws_intersect_periods = reshape(laserON_sws_intersections,2,length(laserON_sws_intersections)/2)'; 
sws_during_laserON_dur = laserON_sws_intersect_periods(:,2)-laserON_sws_intersect_periods(:,1); 
sws_during_laserON_cumdur = sum(sws_during_laserON_dur); 

% total duration of sws during laser OFF periods
laserOFF_sws_intersections = range_intersection(laser_OFFperiods',NREMinclMA_periods_cut'); 
laserOFF_sws_intersect_periods = reshape(laserOFF_sws_intersections,2,length(laserOFF_sws_intersections)/2)'; 
sws_during_laserOFF_dur = laserOFF_sws_intersect_periods(:,2)-laserOFF_sws_intersect_periods(:,1); 
sws_during_laserOFF_cumdur = sum(sws_during_laserOFF_dur); 

oscillation_frq_laserON = laserON_pks_total/sws_during_laserON_cumdur;
oscillation_frq_laserOFF = laserOFF_pks_total/sws_during_laserOFF_cumdur;
