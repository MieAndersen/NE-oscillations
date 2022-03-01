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

   
%% 8) State transitions

% Creating one vector with different behaviors represented by unique
% numbers (1=wake, 4=sws, 9=REM, 15=MA) at frequency 1Hz
boutscore_vector = zeros([1, (Info.HypnoFiles.Duration)+1]);

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

for i=1:length(REM_onset_cut)
    t = REM_onset_cut(i)+1;
    d = REM_duration_cut(i)-1;
    boutscore_vector(t:t+d) = 9; %REM=9
end

for i=1:length(MA_onset_cut)
    t = MA_onset_cut(i)+1;
    d = MA_duration_cut(i)-1;
    boutscore_vector(t:t+d) = 15; %MA=15
end

% Vectors indicate time of transitions in seconds
transition_sws_wake =  find(diff(boutscore_vector)== -3);
transition_wake_sws =  find(diff(boutscore_vector)== 3);
transition_REM_wake =  find(diff(boutscore_vector)== -8);
transition_sws_MA =  find(diff(boutscore_vector)== 11);
transition_REM_sws =  find(diff(boutscore_vector)== -5);
transition_sws_REM =  find(diff(boutscore_vector)== 5);
transition_REM_MA =  find(diff(boutscore_vector)== 6);


%% 9) Plotting all traces and scorings together

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


%% Find LC events based on threshold (Haj's code)
selection = delta3_filt; %the peaks from this trace will be selected and the analysis will be based on these.

% LC Ca2+ event detection. Results are given as indices - not time points
LC_event = DetectAstroEvent(selection,signal_fs, 1.5, true); %2.5 for astrocytes and LC, 1 for neuron (default threshold was 1.5)


%% 10) find LC peaks during sws

period = sws_periods_cut;
period_idx = round(period*signal_fs); %converts time on-/offsets to sampling points in FP dataset

Delta3_stitch = [];
sec_signal_stitch = [];
duration_sws = [];

% Here the sws periods of the full LC trace are extracted and stitched together
for i = 1:length(period_idx)
   period1 = period_idx(i,:);
   trace3 = delta3_filt(period1(1):period1(2)); % for each iteration this will extract a trace excerpt from the full trace for the given period
   time_trace = sec_signal(period1(1):period1(2));
   Delta3_stitch = [Delta3_stitch trace3]; % stitched vector of trace excerpts that lie within the defined period
   sec_signal_stitch = [sec_signal_stitch time_trace]; % a stitch of the time vector so time of peaks can be found
   period_sec = period(i,:);
   sws_time = period_sec(2)-period_sec(1);
   duration_sws = [duration_sws sws_time];
end

% #237: 'MinPeakProminence',2
% #243: 'MinPeakProminence',3.5
minPP = 1.5;   % to use individual thresholds for min preak prominence
%minPP = 4;        % if all recordings use same threshold 

% The stitched trace is used to find peaks occuring during sws
[LC_pks, pklocs, w, p] = findpeaks(Delta3_stitch, sec_signal_stitch, 'MinPeakDistance',5, 'MinPeakProminence',minPP);
    % LC_pks = amplitude of LC peaks, pklocs = time of peaks

% Below I find the time indices of peaks in the full trace
LC_pklocs_idx = zeros(length(pklocs),1);
for i=1:length(pklocs)
        peaktime_idx = find(sec_signal==pklocs(i)); %idx in the full time trace that match peak time
        LC_pklocs_idx(i) = peaktime_idx; % vector of indicex
end

% plotting the events found during sws
% figure 
% findpeaks(Delta3_stitch, sec_signal_stitch, 'MinPeakDistance',5, 'MinPeakProminence',1.5);

figure
    plot(sec_signal(1000:end),delta3_filt(1000:end),'r');
    hold on;
    plot(LC_pklocs_idx./signal_fs, delta3_filt(LC_pklocs_idx),'ko');
    xlabel('time (s)');


%% 11) Triggered LC traces - divided into sws and MA
% LC triggered LC traces

normalize_to_peak = 2.3881; %237 - This is max value of mean sws trace and will be used to normalize (matlab mean: 2.3881; Prism mean: 2.292)
%normalize_to_peak = 6.2017; %243 (matlab mean)

BeforeTime = 15;  %desired number of seconds leading up to LC event (for plotting)
AfterTime = 15;   %desired number of seconds after LC event

BeforeSample = round(BeforeTime * signal_fs); %convert time into # of sampling points
AfterSample = round(AfterTime * signal_fs);

%Extract traces from surrounding LC events
[LC_event_traces,ValidLC_event] = TriggerTraces(delta3_filt,LC_pklocs_idx,BeforeSample,AfterSample);

% time points (in seconds) for LC events 
LC_event_times = sec_signal(ValidLC_event); 

%remove LC events occuring just prior to waking up
min_time2wake = 5; % minimum number of sec from LC event to sws>wake transition

trigger_LC_waking = logical(zeros([1,length(ValidLC_event)])');
for i=1:length(LC_event_times)
    next_transition_wake = transition_sws_wake(find(transition_sws_wake > LC_event_times(i), 1));
    time2wake = next_transition_wake - LC_event_times(i); % time (s)
    if time2wake < min_time2wake % Only LC event less than 5 s prior to sws>MA transition
        trigger_LC_waking(i) = 1;
    end
end

% Which LC events happens less than 5s from sws>MA transition
max_time2transition = 5; % maximum time from LC event to sws>MA transition to count as LC-triggered

trigger_LC_MA = logical(zeros([1,length(ValidLC_event)])');
for i=1:length(LC_event_times)
    next_transition = transition_sws_MA(find(transition_sws_MA > LC_event_times(i), 1));  % this returns the first transition after the LC event
    time2transition = next_transition - LC_event_times(i); % time (s) from LC event to next sws>MA transition
    if time2transition < max_time2transition % Only LC event less than 5 s prior to sws>MA transition
        trigger_LC_MA(i) = 1;
    end
end

% Which of the LC events occurs during sws but not just prior to MAs
[TriggerBool_sws] = IsInInterval(ValidLC_event/signal_fs, sws_periods_cut); %finds LC events occuring during sws
TriggerBool_sws_rmWake = remove_collisions_from_vector(trigger_LC_waking, TriggerBool_sws); % removes LC events from TriggerBool_sws that are classified as trigger_LC_MA
TriggerBool_sws_rmMA = remove_collisions_from_vector(trigger_LC_MA, TriggerBool_sws_rmWake); % removes LC events from TriggerBool_sws that are classified as trigger_LC_MA

% length(ValidLC_event) % total LC events during sws
% sum(trigger_LC_waking) % sws>awake
% sum(TriggerBool_sws_rmMA) %sws
% sum (trigger_LC_MA) %sws>MA
% length(ValidLC_event)-sum(trigger_LC_waking)-sum(TriggerBool_sws_rmMA)-sum (trigger_LC_MA) %undefined event

% select only events that occurs in sws or MA, respectively
sws_LCEventIdx = ValidLC_event(TriggerBool_sws_rmMA);
% sws_LCEventIdx2 = ValidLC_event(trigger_LC_MA==0);    % Why does this give different result?
MA_LCEventIdx = ValidLC_event(trigger_LC_MA);
waking_LCEventIdx = ValidLC_event(trigger_LC_waking);

%Time points (s) of peaks
LC_pktime_sws = sec_signal(sws_LCEventIdx);
LC_pktime_MA = sec_signal(MA_LCEventIdx);
LC_pktime_waking = sec_signal(waking_LCEventIdx);


% group LC traces into sws and MA associated LC events
MA_LC_event = LC_event_traces(:,trigger_LC_MA)/normalize_to_peak;
sws_LC_event = LC_event_traces(:,TriggerBool_sws_rmMA)/normalize_to_peak;
waking_LC_event = LC_event_traces(:,trigger_LC_waking)/normalize_to_peak;
% MA_LC_event = LC_event_traces(:,trigger_LC_MA);
% sws_LC_event = LC_event_traces(:,TriggerBool_sws_rmMA);
% waking_LC_event = LC_event_traces(:,trigger_LC_waking);

    % collect_sws_LC_event = sws_LC_event;
    % collect_MA_LC_event = MA_LC_event;
    % collect_waking_LC_event = waking_LC_event;

    % clearvars -EXCEPT collect_MA_LC_event collect_sws_LC_event collect_waking_LC_event

    % collect_sws_LC_event = horzcat(collect_sws_LC_event,sws_LC_event);
    % collect_MA_LC_event = horzcat(collect_MA_LC_event,MA_LC_event);
    % collect_waking_LC_event = horzcat(collect_waking_LC_event,waking_LC_event);

    % MA_LC_event = collect_MA_LC_event;
    % sws_LC_event = collect_sws_LC_event;
    % waking_LC_event = collect_waking_LC_event;

LC_event_traces_time = ((1:length(MA_LC_event))/signal_fs) - BeforeTime; %time vector in seconds. Time(0) = LC event 

figure
    plot(LC_event_traces_time, mean(MA_LC_event, 2), 'DisplayName','MA')
    hold on
    plot(LC_event_traces_time, mean(sws_LC_event,2), 'DisplayName','sws')
    plot(LC_event_traces_time, mean(waking_LC_event,2), 'DisplayName','waking')
    legend
    xlabel('time (s)');
    ylabel('dF/F');

figure
    plot(LC_event_traces_time, MA_LC_event)
figure
    plot(LC_event_traces_time, sws_LC_event)
figure
    plot(LC_event_traces_time, waking_LC_event)

prism_downsampling = 10;

% mean trace
mean_LC_MA = mean(MA_LC_event, 2);
mean_LC_MA_ds = downsample(mean_LC_MA, prism_downsampling);
mean_LC_sws = mean(sws_LC_event,2);
mean_LC_sws_ds = downsample(mean_LC_sws, prism_downsampling);
mean_LC_waking = mean(waking_LC_event, 2);
mean_LC_waking_ds = downsample(mean_LC_waking, prism_downsampling);


% SEM of traces
SEM_LC_MA = downsample(std(MA_LC_event, [], 2)./ sqrt(size(MA_LC_event,2)),prism_downsampling);
SEM_LC_sws = downsample(std(sws_LC_event, [], 2)./ sqrt(size(sws_LC_event,2)),prism_downsampling);
SEM_LC_waking = downsample(std(waking_LC_event, [], 2)./ sqrt(size(waking_LC_event,2)),prism_downsampling);

% timeline for prism - starting at 0
prism_timeline = downsample(LC_event_traces_time, prism_downsampling)';
prism_timeline_zero = (prism_timeline-prism_timeline(1));

% input vector for n (for prism)
number_vector_MA = zeros(1,length(prism_timeline))'+size(MA_LC_event,2);
number_vector_sws = zeros(1,length(prism_timeline))'+size(sws_LC_event,2);
number_vector_waking = zeros(1,length(prism_timeline))'+size(waking_LC_event,2);


% Amplitude of peaks
peakVal_LC_MA = LC_pks(trigger_LC_MA)'/normalize_to_peak;
peakVal_LC_sws = LC_pks(TriggerBool_sws_rmMA)'/normalize_to_peak;
peakVal_LC_waking = LC_pks(trigger_LC_waking)'/normalize_to_peak;
% peakVal_LC_MA = LC_pks(trigger_LC_MA;
% peakVal_LC_sws = LC_pks(TriggerBool_sws_rmMA)';
% peakVal_LC_waking = LC_pks(trigger_LC_waking)';

% peak prominence 
peakprom_LC_MA = p(trigger_LC_MA)'/normalize_to_peak;
peakprom_LC_sws = p(TriggerBool_sws_rmMA)'/normalize_to_peak;
peakprom_LC_waking = p(trigger_LC_waking)'/normalize_to_peak;
% peakprom_LC_MA = p(trigger_LC_MA)';
% peakprom_LC_sws = p(TriggerBool_sws_rmMA)';
% peakprom_LC_waking = p(trigger_LC_waking)';

% baseline before/after LC event
baseline_LC_MA = peakVal_LC_MA-peakprom_LC_MA;
baseline_LC_sws = peakVal_LC_sws-peakprom_LC_sws;
baseline_LC_waking = peakVal_LC_waking-peakprom_LC_waking;


% downsample individual traces
MA_LC_event_ds = downsample(MA_LC_event, prism_downsampling);
sws_LC_event_ds = downsample(sws_LC_event, prism_downsampling);
waking_LC_event_ds = downsample(waking_LC_event, prism_downsampling);

% Find the time from LC event to MA onset
    % LC event times (s) of events that lead to MAs
    LC_event_times_MA = LC_event_times(trigger_LC_MA);
    LC_event_times_waking = LC_event_times(trigger_LC_waking); 
    
    % Find the onset of MA for each LC triggered MA
    transition_sws_MA_LC = zeros(1,length(LC_event_times_MA))'; 
    for i=1:length(LC_event_times_MA)
        next_transition = transition_sws_MA(find(transition_sws_MA > LC_event_times_MA(i), 1));  % this returns the first transition after the LC event
        transition_sws_MA_LC(i) = next_transition; % onset times (s) of MAs triggered by LC event
    end
    
    % Find the onset of waking for each LC triggered waking
    transition_waking_LC = zeros(1,length(LC_event_times_waking))'; 
    for i=1:length(LC_event_times_waking)
        next_transition = transition_sws_wake(find(transition_sws_wake > LC_event_times_waking(i), 1));  % this returns the first transition after the LC event
        transition_waking_LC(i) = next_transition; % onset times (s) of MAs triggered by LC event
    end    
    
    
    % time from LC event to associated MA/wake onset
    LCevent2MA_time = transition_sws_MA_LC - LC_event_times_MA';
    LCevent2waking_time = transition_waking_LC - LC_event_times_waking'; 
    
%     collect_LCevent2MA_time = LCevent2MA_time;
%     collect_LCevent2waking_time = LCevent2waking_time;

%     clearvars -EXCEPT collect_LCevent2MA_time collect_LCevent2waking_time

%     collect_LCevent2MA_time = vertcat(collect_LCevent2MA_time,LCevent2MA_time);
%     collect_LCevent2waking_time = vertcat(collect_LCevent2waking_time,LCevent2waking_time);
% 
%     LCevent2MA_time = collect_LCevent2MA_time;
%     LCevent2waking_time = collect_LCevent2waking_time;

    
    % mean and SEM of time from LC to MA
    meanLCevent2MA_time = mean(LCevent2MA_time);
    SEMLCevent2MA_time = std(LCevent2MA_time)./ sqrt(length(LCevent2MA_time));
    
    meanLCevent2waking_time = mean(LCevent2waking_time);
    SEMLCevent2waking_time = std(LCevent2waking_time)./ sqrt(length(LCevent2waking_time));
    
    
% Optimizing minPP
    % length(LC_event_times); % #LC events
    % length(LC_event_times_MA); % #LC-MA
    % length(sws_LCEventIdx); % #LC-sws
    % length(LC_event_times)-length(LC_event_times_MA)-length(sws_LCEventIdx); % #undef
    % length(MA_duration); % #MA

% AUC
% to remove negative values from the AUC calculation, a constant (add_auc) will be added to all dF/F values, 
lowest_dF = min([min(min(MA_LC_event)), min(min(sws_LC_event)), min(min(waking_LC_event))]);
if lowest_dF < 0
    add_to_dF = abs(lowest_dF);
else
    add_to_dF = 0;
end

% AUC time
auc_start = 2; %AUC will be calculated from this time (s) relative to time=0 
auc_end = 10;   %and until this time (s) relative to time=0
auc_start_idx = find(LC_event_traces_time>auc_start,1);
auc_end_idx = find(LC_event_traces_time>auc_end,1);

% MA AUC
AUC_preLC_event_MA = zeros(size(MA_LC_event,2),1);
for i=1:size(MA_LC_event,2)
    auc_i = trapz(LC_event_traces_time(auc_start_idx:auc_end_idx), MA_LC_event(auc_start_idx:auc_end_idx,i)'+add_to_dF);
    AUC_preLC_event_MA(i,1) = auc_i;
end

%sws AUC
AUC_preLC_event_sws = zeros(size(sws_LC_event,2),1);
for i=1:size(sws_LC_event,2)
    auc_i = trapz(LC_event_traces_time(auc_start_idx:auc_end_idx), sws_LC_event(auc_start_idx:auc_end_idx,i)'+add_to_dF);
    AUC_preLC_event_sws(i,1) = auc_i;
end

% waking AUC
AUC_preLC_event_waking = zeros(size(waking_LC_event,2),1);
for i=1:size(waking_LC_event,2)
    auc_i = trapz(LC_event_traces_time(auc_start_idx:auc_end_idx), waking_LC_event(auc_start_idx:auc_end_idx,i)'+add_to_dF);
    AUC_preLC_event_waking(i,1) = auc_i;
end


%% 13) LC Triggered NE traces

% LC triggered NE traces

%normalize_to_NEpeak = 0.7068 ; %237 - This is amplitude of mean sws trace and will be used to normalize
normalize_to_NEpeak = 1.3894 ; %243

BeforeTime = 60;  %desired number of seconds leading up to LC event (for plotting)
AfterTime = 120;   %desired number of seconds after LC event

BeforeSample = round(BeforeTime * signal_fs); %convert time into # of sampling points
AfterSample = round(AfterTime * signal_fs);

%Extract NE traces before/after LC event during sws
[LC_event_NEtraces,ValidNE_event] = TriggerTraces(delta2_filt,LC_pklocs_idx,BeforeSample,AfterSample);

% divide traces into sws and MA related 
MA_NE_LC_event = LC_event_NEtraces(:,trigger_LC_MA)/normalize_to_NEpeak;
sws_NE_LC_event = LC_event_NEtraces(:,TriggerBool_sws_rmMA(1:end-1))/normalize_to_NEpeak;
waking_NE_LC_event = LC_event_NEtraces(:,trigger_LC_waking)/normalize_to_NEpeak;
%Non-normalized:
% MA_NE_LC_event = LC_event_NEtraces(:,trigger_LC_MA); 
% sws_NE_LC_event = LC_event_NEtraces(:,TriggerBool_sws_rmMA);
% waking_NE_LC_event = LC_event_NEtraces(:,trigger_LC_waking);


%     collect_sws_NE_event = sws_NE_LC_event;
%     collect_MA_NE_event = MA_NE_LC_event;
%     collect_waking_NE_event = waking_NE_LC_event;

%     clearvars -EXCEPT collect_MA_NE_event collect_sws_NE_event collect_waking_NE_event
% 
%     collect_sws_NE_event = horzcat(collect_sws_NE_event,sws_NE_LC_event);
%     collect_MA_NE_event = horzcat(collect_MA_NE_event,MA_NE_LC_event);
%     collect_waking_NE_event = horzcat(collect_waking_NE_event,waking_NE_LC_event);
% 
%     MA_NE_LC_event = collect_MA_NE_event;
%     sws_NE_LC_event = collect_sws_NE_event;
%     waking_NE_LC_event = collect_waking_NE_event;


NE_traces_time = ((1:length(MA_NE_LC_event))/signal_fs) - BeforeTime; %time vector in seconds. Time(0) = LC event 

%plot mean NE traces triggered by LC
figure
    plot(NE_traces_time, mean(MA_NE_LC_event, 2), 'DisplayName','MA')
    hold on
    plot(NE_traces_time, mean(sws_NE_LC_event, 2), 'DisplayName','sws')
    plot(NE_traces_time, mean(waking_NE_LC_event, 2), 'DisplayName','waking')
    legend
    xlabel('time (s)');
    ylabel('dF/F');

% plot all NE traces triggered by LC
figure
plot(NE_traces_time, MA_NE_LC_event)
figure
plot(NE_traces_time, sws_NE_LC_event)
figure
plot(NE_traces_time, waking_NE_LC_event)

prism_downsampling = 60;

% mean trace
mean_NE_MA = mean(MA_NE_LC_event, 2);
mean_NE_MA_ds = downsample(mean_NE_MA, prism_downsampling);
mean_NE_sws = mean(sws_NE_LC_event, 2);
mean_NE_sws_ds = downsample(mean_NE_sws, prism_downsampling);
mean_NE_waking = mean(waking_NE_LC_event, 2);
mean_NE_waking_ds = downsample(mean_NE_waking, prism_downsampling);

% SEM of traces
SEM_NE_MA = downsample(std(MA_NE_LC_event, [], 2)./ sqrt(size(MA_NE_LC_event,2)),prism_downsampling);
SEM_NE_sws = downsample(std(sws_NE_LC_event, [], 2)./ sqrt(size(sws_NE_LC_event,2)),prism_downsampling);
SEM_NE_waking = downsample(std(waking_NE_LC_event, [], 2)./ sqrt(size(waking_NE_LC_event,2)),prism_downsampling);

% input vector for n (for prism)
number_vector_MA = zeros(1,length(NE_traces_time))'+size(MA_NE_LC_event,2);
number_vector_sws = zeros(1,length(NE_traces_time))'+size(sws_NE_LC_event,2);
number_vector_waking = zeros(1,length(NE_traces_time))'+size(waking_NE_LC_event,2);

% timeline for prism - starting at 0
prism_timeline = downsample(NE_traces_time, prism_downsampling)';
prism_timeline_zero = (prism_timeline-prism_timeline(1));

%     temp_sws_base = mean(mean_NE_sws(find(NE_traces_time>-6,1):find(NE_traces_time>-1,1)));
%     temp_sws_max = max(mean_NE_sws);
%     temp_sws_ampl = temp_sws_max-temp_sws_base;

% Trace characteristics of NE traces
kinetics_MA = zeros(size(MA_NE_LC_event,2),5);
for i=1:size(MA_NE_LC_event,2)
    [ kinetics_MA(i,1),kinetics_MA(i,2), kinetics_MA(i,3), kinetics_MA(i,4), kinetics_MA(i,5)] = Kinetics_Haj_noPlot( MA_NE_LC_event(:,i), MeanFilterOrder, signal_fs, 15,BeforeTime );
end    

kinetics_sws = zeros(size(sws_NE_LC_event,2),5);
for i=1:size(sws_NE_LC_event,2)
    [ kinetics_sws(i,1),kinetics_sws(i,2), kinetics_sws(i,3), kinetics_sws(i,4), kinetics_sws(i,5)] = Kinetics_Haj_noPlot( sws_NE_LC_event(:,i), MeanFilterOrder, signal_fs, 15,BeforeTime );
% for i=11:20 %check some om of peaks
%   [ kinetics_sws(i,1),kinetics_sws(i,2), kinetics_sws(i,3), kinetics_sws(i,4), kinetics_sws(i,5)] = Kinetics_Haj( sws_NE_LC_event(:,i), MeanFilterOrder, signal_fs, 15,BeforeTime );
end   

kinetics_waking = zeros(size(waking_NE_LC_event,2),5);
for i=1:size(waking_NE_LC_event,2)
    [ kinetics_waking(i,1),kinetics_waking(i,2), kinetics_waking(i,3), kinetics_waking(i,4), kinetics_waking(i,5)] = Kinetics_Haj_noPlot( waking_NE_LC_event(:,i), MeanFilterOrder, signal_fs, 15,BeforeTime );
end   

% Characterististics of mean traces
%[ kinetics_MA(1,1),kinetics_MA(1,2), kinetics_MA(1,3), kinetics_MA(1,4), kinetics_MA(1,5)] = Kinetics_Haj( mean_NE_MA, MeanFilterOrder, signal_fs, 15,BeforeTime);
%[ kinetics_sws(1,1),kinetics_sws(1,2), kinetics_sws(1,3), kinetics_sws(1,4), kinetics_sws(1,5)] = Kinetics_Haj( mean_NE_sws, MeanFilterOrder, signal_fs, 15,BeforeTime );
%[ kinetics_waking(1,1),kinetics_waking(1,2), kinetics_waking(1,3), kinetics_waking(1,4), kinetics_waking(1,5)] = Kinetics_Haj( mean_NE_waking, MeanFilterOrder, signal_fs, 15,BeforeTime );
% outputs: PeakVal,PeakTimeInSmp, Baseline, PeakMin, PeakTime_min, Base_end, RaiseTime20_80, DecayTime80_20, trace_out

% Amplitude of peaks
peakVal_NE_MA = kinetics_MA(:,1);
peakVal_NE_sws = kinetics_sws(:,1);
peakVal_NE_waking = kinetics_waking(:,1);

% baseline before/after LC event
baseline_NE_MA = kinetics_MA(:,3);
baseline_NE_sws = kinetics_sws(:,3);
baseline_NE_waking = kinetics_waking(:,3);

% peak prominence 
peakprom_NE_MA = peakVal_NE_MA-baseline_NE_MA;
peakprom_NE_sws = peakVal_NE_sws-baseline_NE_sws;
peakprom_NE_waking = peakVal_NE_waking-baseline_NE_waking;


% baseline, peak, and amplitude
% NE baseline pre REM offset
NE_bl_start_time = 14; %sec prior to transition
NE_bl_end_time = 1; %sec prior to transition
NE_bl_start_idx = abs(round((NE_bl_start_time-BeforeTime)*signal_fs)); % index of timepoint 5 s prior to peak time
NE_bl_end_idx = abs(round((NE_bl_end_time-BeforeTime)*signal_fs)); % index of timepoint 3 s prior to peak time

%NE baseline prior to LC event
NEbase_preLCevent_sws = zeros(size(sws_NE_LC_event,2),1);
for i=1:size(sws_NE_LC_event,2)
    NEbase_preLC_sws_i = mean(sws_NE_LC_event(NE_bl_start_idx:NE_bl_end_idx,i));
    NEbase_preLCevent_sws(i,1) = NEbase_preLC_sws_i;
end

NEbase_preLCevent_MA = zeros(size(MA_NE_LC_event,2),1);
for i=1:size(MA_NE_LC_event,2)
    NEbase_preLC_MA_i = mean(MA_NE_LC_event(NE_bl_start_idx:NE_bl_end_idx,i));
    NEbase_preLCevent_MA(i,1) = NEbase_preLC_MA_i;
end

NEbase_preLCevent_waking = zeros(size(waking_NE_LC_event,2),1);
for i=1:size(waking_NE_LC_event,2)
    NEbase_preLC_waking_i = mean(waking_NE_LC_event(NE_bl_start_idx:NE_bl_end_idx,i));
    NEbase_preLCevent_waking(i,1) = NEbase_preLC_waking_i;
end

%peak time for the mean trace is used to decide where to look for peak
%values in the individual traces
[~,NEpeak_idx_sws] = max(mean_NE_sws);
[~,NEpeak_idx_MA] = max(mean_NE_MA);
[~,NEpeak_idx_waking] = max(mean_NE_waking);
sec_after_meanPeak = round(10*signal_fs);

% NE peak LC event
NEpeak_postLC_sws = zeros(size(sws_NE_LC_event,2),1);
for i=1:size(sws_NE_LC_event,2)
    %NEpeak_postLC_sws_i = max(sws_NE_LC_event(BeforeSample:NEpeak_idx_sws+sec_after_meanPeak,i)); %finds the maximum value in the NE trace after REM offset and [AfterTime] seconds after.
    NEpeak_postLC_sws_i = max(sws_NE_LC_event(BeforeSample:end,i));
    NEpeak_postLC_sws(i,1) = NEpeak_postLC_sws_i;
end

NEpeak_postLC_MA = zeros(size(MA_NE_LC_event,2),1);
for i=1:size(MA_NE_LC_event,2)
    %NEpeak_postLC_MA_i = max(MA_NE_LC_event(BeforeSample:NEpeak_idx_MA+sec_after_meanPeak,i)); %finds the maximum value in the NE trace after REM offset and [AfterTime] seconds after.
    NEpeak_postLC_MA_i = max(MA_NE_LC_event(BeforeSample:end,i));
    NEpeak_postLC_MA(i,1) = NEpeak_postLC_MA_i;
end

NEpeak_postLC_waking = zeros(size(waking_NE_LC_event,2),1);
for i=1:size(waking_NE_LC_event,2)
    %NEpeak_postLC_waking_i = max(waking_NE_LC_event(BeforeSample:NEpeak_idx_waking+sec_after_meanPeak,i)); %finds the maximum value in the NE trace after REM offset and [AfterTime] seconds after.
    NEpeak_postLC_waking_i = max(waking_NE_LC_event(BeforeSample:end,i));
    NEpeak_postLC_waking(i,1) = NEpeak_postLC_waking_i;
end

% NE amplitude
NEampl_LCevent_sws = NEpeak_postLC_sws-NEbase_preLCevent_sws;
NEampl_LCevent_MA = NEpeak_postLC_MA-NEbase_preLCevent_MA;
NEampl_LCevent_waking = NEpeak_postLC_waking-NEbase_preLCevent_waking;


%baseline, peak, and amplitude of mean trace
NEbase_sws_mean = mean(mean_NE_sws(NE_bl_start_idx:NE_bl_end_idx));
NEpeak_sws_mean = max(mean_NE_sws(BeforeSample:end));
NEampl_ses_mean = NEpeak_sws_mean-NEbase_sws_mean; % this value is used to normalize for each animal (normalize_to_NEpeak)

%downsampling indvdl traces for plotting
MA_NE_LC_event_ds = downsample(MA_NE_LC_event, prism_downsampling);
sws_NE_LC_event_ds = downsample(sws_NE_LC_event, prism_downsampling);
waking_NE_LC_event_ds = downsample(waking_NE_LC_event, prism_downsampling);


% to remove negative values from the AUC calculation, a constant (add_auc) will be added to all dF/F values, 
lowest_dF = min([min(min(MA_NE_LC_event)), min(min(sws_NE_LC_event)), min(min(waking_NE_LC_event))]);
if lowest_dF < 0
    add_to_dF = abs(lowest_dF);
else
    add_to_dF = 0;
end

% AUC time
auc_start = -30; %AUC will be calculated from this time (s)
auc_end = 30;   %and until this time (s)
auc_start_idx = find(NE_traces_time>auc_start,1);
auc_end_idx = find(NE_traces_time>auc_end,1);

% MA AUC
AUC_NE_preLC_event_MA = zeros(size(MA_NE_LC_event,2),1);
for i=1:size(MA_NE_LC_event,2)
    auc_i = trapz(NE_traces_time(auc_start_idx:auc_end_idx), MA_NE_LC_event(auc_start_idx:auc_end_idx,i)'+add_to_dF);
    AUC_NE_preLC_event_MA(i,1) = auc_i;
end

%sws AUC
AUC_NE_preLC_event_sws = zeros(size(sws_NE_LC_event,2),1);
for i=1:size(sws_NE_LC_event,2)
    auc_i = trapz(NE_traces_time(auc_start_idx:auc_end_idx), sws_NE_LC_event(auc_start_idx:auc_end_idx,i)'+add_to_dF);
    AUC_NE_preLC_event_sws(i,1) = auc_i;
end

% waking AUC
AUC_NE_preLC_event_waking = zeros(size(waking_NE_LC_event,2),1);
for i=1:size(waking_NE_LC_event,2)
    auc_i = trapz(NE_traces_time(auc_start_idx:auc_end_idx), waking_NE_LC_event(auc_start_idx:auc_end_idx,i)'+add_to_dF);
    AUC_NE_preLC_event_waking(i,1) = auc_i;
end

% figure
% plot(NE_traces_time(auc_start_idx:auc_end_idx), waking_NE_LC_event(auc_start_idx:auc_end_idx,i)'+add_to_dF)

%% 14) REM 

% MUST ONLY BE RUN ONCE!!!

%manual removal of questionable REM bouts

%{
keep_REMbouts = logical(ones(length(REM_periods_cut),1));
rm_REMbouts = [2 3 4 11 12];                                          % <<---- list bout# (row numbers) to be removed (2, 3, 4, 11, 12 for #237, and 1 for #243) 
keep_REMbouts(rm_REMbouts)=0;
REM_periods_cut = REM_periods_cut(keep_REMbouts,:);
REM_offset_cut = REM_offset_cut(keep_REMbouts,:);
REM_onset_cut = REM_onset_cut(keep_REMbouts,:);
REM_duration_cut = REM_duration_cut(keep_REMbouts,:);
%}

%Manual scoring of REM onsets/offsets based on NE trace - but only for bouts recognized in EEG score
EEG_FP_triple_237_manualREMscore = [2291 2365; 2516 2632; 2877 3024; 4603 4780; 6642 6870; 7461 7567; 8078 8244; 9746 9965; 10828 10918; 11360 11436; 11592 11649; 11743 11840];
EEG_FP_triple_243_manualREMscore = [5585 5780; 7812 7969; 8832 8959; 10441 10665];


%% 15) REM drop

% REM drop vs duration
period = REM_periods_cut;
period_fs = period*signal_fs;

% Delta1_all = [];
% Delta2_all = [];
% Delta3_all = [];
REM_duration2 = zeros(1,length(period_fs));
REM_drop = zeros(1,length(period_fs));
REM_start = zeros(1,length(period_fs));

for i = 1:length(period_fs)
   period1 = period_fs(i,:);
   period_sec = period(i,:);
%    trace1 = delta1_filt(period1(1):period1(2));
   trace2 = delta2_filt(period1(1):period1(2));
%    trace3 = delta3_filt(period1(1):period1(2));
   REM_time = period_sec(2)-period_sec(1);
   %REM_base = trace2(end) - trace2(1);
   REM_base = min(trace2) - max(trace2);
%    Delta1_all = [Delta1_all trace1];
%    Delta2_all = [Delta2_all trace2];
%    Delta3_all = [Delta3_all trace3];
   REM_duration2(i) = REM_time;
   REM_drop(i) = REM_base;  % magnitude of the drop drom basline
   REM_start(i) = mean(trace2(1:100)); 
end

%plot REM drop vs duration
figure
    scatter(REM_duration2,REM_drop)
    xlabel('REM duration (s)');
    ylabel('REM drop (dF/F)');

mean_REM_duration = mean(REM_duration2);
mean_REM_drop = mean(REM_drop);
mean_REM_start = mean(REM_start);
REM_occurences = length(period_fs)/(length(delta1_filt)/signal_fs/60/60);

% Manually defined REM drop on-/offsets
EEG_FP_237_REM_NE = ([2291 2365; 2516 2632; 2873 3024; 4603 4780; 6642 6870; 7461 7567; 8078 8244; 9742 9965; 10828 10918; 11360 11436; 11592 11649; 11743 11840]); %the 100 sec were removed from trace to align with FP 4s were removed to aign with socre - consider removing here as well.
EEG_FP_243_REM_NE = ([5585 5780; 7812 7969; 8832 8959; 10441 10665]); %the 58 sec were removed from trace to align with FP. 4s were added to align with score - consider removing here as well.
EEG_FP_329_REM_NE =([6084.62 6281.53; 8961.39 9114.16]);
EEG_FP_323_REM_NE = ([3967.44 4104.24; 5817.72 5932.88; 6798.92 7092.46]);
EEG_FP_319_REM_NE = ([3206.2 3416.6; 4483.6 4710.2; 7363.5 7446.5; 7578.7 7690.8]);
EEG_FP_331_REM_NE = [3532.16 3753.44; 4776.28 4881.31; 6602.73 6804.21; 7374.13 7466.64; 7920.31 8022.65];
EEG_FP_307_REM_NE =[3852.46 4100.68; 5752.66 5910.98; 7820.73 8020.58; 9452.03 9686.74; 11810.69 12064.76];

REM_drop_periods = EEG_FP_329_REM_NE;

% Find decay time for REM drop
REMdrop_tau = zeros(length(REM_onset_cut),1);
REMdrop_halflife = zeros(length(REM_onset_cut),1);

remove_tail = 0; % # of seconds removed from REM tail for fitting purpose
add_tail = 5; % add tail to REM bout for fitting purposes
remove_start = 0; % remove first # of seconds of REM bout for fitting purposes
REM_smooth = 1000;

% Decay for individual REM bouts
for i= 2 %1:length(REM_onset_cut)
    REM_trace_start = round((REM_drop_periods(i,1)+remove_start)*signal_fs); % index of REM onset in full trace
    REM_trace_end = round((REM_drop_periods(i,2)+add_tail)*signal_fs); % index of REM offset in full trace
    %REM_trace_start = round(REM_periods_cut(i,1)*signal_fs); % index of REM onset in full trace
    %REM_trace_end = round((REM_periods_cut(i,2)-remove_tail)*signal_fs); % index of REM offset in full trace
    REM_trace = smooth(delta2_filt(REM_trace_start:REM_trace_end),REM_smooth); % trace extract
    REM_trace_time = (1:length(REM_trace))/signal_fs;
    fit_type = fittype('a+b*exp(-c*x)');
    StartPoint = [min(REM_trace(:)),0,0]; % because values go way below zero we want to 
    f = fit(REM_trace_time(:), REM_trace(:), fit_type,'StartPoint',StartPoint); % matlab will optimize the fit based on the startpoint 
    figure
        plot(f,REM_trace_time', REM_trace)
    t_half = log(2)/f.c; %halflife. (c=lambda is the decay constant)
    time_constant = 1/(f.c);  %tau (time constant) = 1/lambda
    REMdrop_tau(i) = time_constant;
    REMdrop_halflife(i) = t_half;
end


%plot REM drop vs duration
figure
    scatter(REMdrop_tau,REM_drop)
    xlabel('tau');
    ylabel('REM drop (dF/F)');
figure
    scatter(REMdrop_halflife,REM_drop)
    xlabel('halflife');
    ylabel('REM drop (dF/F)');


% Decay for mean of REM bouts (initial 30 s)

before = 20; % from onset of NE drop
after = 0; % s after NE drop

% Manually exclude onsets that are off
    exclude_N2REM = [NaN]; % <<-------------- PUT IDX OF ONSETS TO BE REMOVED (put NaN if none should be removed)
    REMdrop_onset_xcld = REM_periods_cut(:,1);
    REM_period_xcld = REM_periods_cut;
    REMdrop_onset_xcld(exclude_N2REM) = []; %skip if none are removed
    REM_period_xcld(exclude_N2REM,:) = [];
    
    % remove drops shorter than trace extracts
    REMdrop_onset_subset = REMdrop_onset_xcld(REM_period_xcld(:,2)-REMdrop_onset_xcld>after); % remove drops shorter than 20 s

    % Trace extracts of intial part of REM drop
    REM_period_epocs2 = [];
    REM_period_epocs3 = [];
    for i = 1:length(REMdrop_onset_subset) % using REM onset
         time = REMdrop_onset_subset(i);  
         signal2_epoc = delta2_filt(round((time- before)*signal_fs:(time + after)*signal_fs));
         signal3_epoc = delta3_filt(round((time- before)*signal_fs:(time + after)*signal_fs));
         REM_period_epocs2(:,i) = signal2_epoc;
         REM_period_epocs3(:,i) = signal3_epoc;
    end

    fs_signal = 1:1:length(signal2_epoc);
    sec_signal = (fs_signal/signal_fs)-before;

    % mean traces from initial drop
    figure
    a = subplot(2,1,1);
        plot(downsample(sec_signal,20), downsample(mean(REM_period_epocs2,2), 20));    
        title('REM onset - NE')
    b = subplot(2,1,2);
        plot(downsample(sec_signal,20), downsample(mean(REM_period_epocs3,2), 20));    
        title('REM onset - LC')
    linkaxes([a,b],'x');

% Mean REM drop trace
mean_REM_period_epocs2 = mean(REM_period_epocs2,2); % all REM pooled
ds_mean_REM_period_epocs2 = downsample(mean_REM_period_epocs2,20);
ds_sec_signal = downsample(sec_signal,20)';


% Fit exponential model
%{
fit_trace = mean_REM_period_epocs2;

fit_type = fittype('a+b*exp(-c*x)');
StartPoint = [min(fit_trace(:)),0,0]; % because values go way below zero we want to 
%f = fit(sec_signal(:), fit_trace(:),fit_type,'StartPoint',StartPoint); % matlab will optimize the fit based on the startpoint 
f = fit(sec_signal(:), fit_trace(:),fit_type); % 
figure
    plot(f,sec_signal', fit_trace)
t_half = log(2)/f.b; %halflife. (b=lambda is the decay constant)
t_constant = 1/(f.b);  %tau (time constant) = 1/b
%}
  
    
%% 16) Episode following REM

    %{
    OBS!
    To run this section, awake bouts must NOT be subdivided into wake vs
    MA. Skip section 6 and 8 when running the code - and exclude MA related
    lines in section 7. If, this has already been done, run section 5 again
    to reload original sleep score.
    %}

% Find which wake bouts that follow REM bouts
wake_follow_REM_idx = zeros(length(REM_offset_cut),1);

for i=1:length(REM_offset_cut)
    following_wake_idx = find(REM_offset_cut(i)<=wake_onset_cut,1); % this gives the index of the coninsiding or first following wake bout
    
    % consider including if statement so wake bouts with onset more than 5s after REM offset are put as NaN 
    
    wake_follow_REM_idx(i) = following_wake_idx;
end

wake_follow_REM_duration = wake_duration_cut(wake_follow_REM_idx);

figure
    scatter(wake_follow_REM_duration, REM_drop')
    xlabel('following wake bout duration (s)');
    ylabel('REM drop (dF/F)');
figure
    scatter(wake_follow_REM_duration, REM_duration_cut)
    xlabel('following wake bout duration (s)');
    ylabel('REM duration (s)');

figure
    scatter(wake_follow_REM_duration, REMdrop_tau)
    xlabel('following wake bout duration (s)');
    ylabel('tau');
figure
    scatter(wake_follow_REM_duration, REMdrop_halflife)
    xlabel('following wake bout duration (s)');
    ylabel('halflife');


%type of episode


%% 17) REM onset traces

BeforeTime = 60;  %desired number of seconds leading up to REM onset (for plotting)
AfterTime = 60;   %desired number of seconds after REM onset
    %for NE use 10s before, 30s after to get drop
    %for LC: depends
prism_downsampling = 30;

BeforeSample = round(BeforeTime * signal_fs); %convert time into # of sampling points
AfterSample = round(AfterTime * signal_fs);

% Below I find the time indices of REM onsets - to index into
REM_onset_idx = round(REM_onset_cut*signal_fs);

%Extract NE and LC traces before/after REM onset
[REM_onset_LCtraces,~] = TriggerTraces(delta3_filt,REM_onset_idx,BeforeSample,AfterSample);
[REM_onset_NEtraces,ValidNE_event] = TriggerTraces(delta2_filt,REM_onset_idx,BeforeSample,AfterSample);

%downsampled traces
REM_onset_LCtraces_ds = downsample(REM_onset_LCtraces, prism_downsampling);
REM_onset_NEtraces_ds = downsample(REM_onset_NEtraces, prism_downsampling);

REMon_traces_time = ((1:length(REM_onset_LCtraces))/signal_fs) - BeforeTime; %time vector in seconds. Time(0) = LC event 

%plot mean LC and NE traces at REM onset
figure
    plot(REMon_traces_time, mean(REM_onset_LCtraces, 2), 'DisplayName','LC')
    hold on
    plot(REMon_traces_time, mean(REM_onset_NEtraces,2), 'DisplayName','NE')
    legend
    xlabel('time (s)');
    ylabel('dF/F');

% plot all NE traces triggered by LC
figure
plot(REMon_traces_time, REM_onset_NEtraces)
figure
plot(REMon_traces_time, REM_onset_LCtraces)

% mean trace
mean_LC_REMonset = mean(REM_onset_LCtraces, 2);
mean_LC_REMonset_ds = downsample(mean_LC_REMonset, prism_downsampling);
mean_NE_REMonset = mean(REM_onset_NEtraces,2);
mean_NE_REMonset_ds = downsample(mean_NE_REMonset, prism_downsampling);

% SEM of traces
SEM_LC_REMon = downsample(std(REM_onset_LCtraces, [], 2)./ sqrt(size(REM_onset_LCtraces,2)),prism_downsampling);
SEM_NE_REMon = downsample(std(REM_onset_NEtraces, [], 2)./ sqrt(size(REM_onset_NEtraces,2)),prism_downsampling);

% input vector for n (for prism)
number_vector_LC_REMon = zeros(1,length(SEM_LC_REMon))'+size(REM_onset_LCtraces,2);
number_vector_NE_REMon = zeros(1,length(SEM_NE_REMon))'+size(REM_onset_NEtraces,2);

% timeline for prism - starting at 0
prism_timeline = downsample(REMon_traces_time, prism_downsampling)';
prism_timeline_zero = (prism_timeline-prism_timeline(1));

% NE decay time at REM onset (of mean trace)
fit_type = fittype('a+b*exp(-c*x)');
StartPoint = [min(mean_NE_REMonset(:)),0,0];
f_REMon = fit(REMon_traces_time', mean_NE_REMonset,fit_type, 'StartPoint',StartPoint); 

time_constant = 1/(f_REMon.b); %tau is 1/b
decay_pos = abs(time_constant);

figure
plot(f_REMon,REMon_traces_time', mean_NE_REMonset)


% LC baseline pre/post REM onset
preREMon_LC_bl_time = 60; %sec prior to transition
postREMon_LC_bl_time = 60; %sec after transition
BeforeSample = round(preREMon_LC_bl_time * signal_fs); %convert time into # of sampling points
AfterSample = round(postREMon_LC_bl_time * signal_fs);

% in this loop 
preREMon_LCtraces = NaN(round(preREMon_LC_bl_time*signal_fs)+1,length(REM_onset_cut));
postREMon_LCtraces = NaN(round(postREMon_LC_bl_time*signal_fs)+1,length(REM_onset_cut));
for i=1:length(REM_onset_cut)
    if REM_duration_cut(i) < preREMon_LC_bl_time % for REM bouts shorter than 1 min baseline before/after are shorter
        preREMon_LC_bl_i = REM_duration_cut(i);
        postREMon_LC_bl_i = REM_duration_cut(i);
        BeforeSample_i = round(preREMon_LC_bl_i * signal_fs); %convert time into # of sampling points
        AfterSample_i = round(postREMon_LC_bl_i * signal_fs);
        preREMon_LCtraces(BeforeSample-BeforeSample_i+1:end,i) = TriggerTraces(delta3_filt, REM_onset_idx(i), BeforeSample_i, 0);
        postREMon_LCtraces(AfterSample-AfterSample_i+1:end,i) = TriggerTraces(delta3_filt, REM_onset_idx(i), 0, AfterSample_i);
    else
        preREMon_LCtraces(:,i) = TriggerTraces(delta3_filt, REM_onset_idx(i), BeforeSample, 0);
        postREMon_LCtraces(:,i) = TriggerTraces(delta3_filt, REM_onset_idx(i), 0, AfterSample);
    end
end

% preREMon_LCtraces = TriggerTraces(delta3_filt, REM_onset_idx, BeforeSample, 0);
% postREMon_LCtraces = TriggerTraces(delta3_filt, REM_onset_idx, 0, AfterSample);

%mean LC baseline
preREMon_LCbl = nanmean(preREMon_LCtraces);
postREMon_LCbl = nanmean(postREMon_LCtraces);

%SEM LC basline
preREMon_LCbl_SEM = nanstd(preREMon_LCtraces', [], 2)./ sqrt(size(preREMon_LCtraces',2));
postREMon_LCbl_SEM = nanstd(postREMon_LCtraces', [], 2)./ sqrt(size(postREMon_LCtraces',2));

number_vector_LC_REMon = zeros(1,length(preREMon_LCbl))'+size(postREMon_LCtraces',2);


%% 18) REM offset traces

BeforeTime = 60;  %desired number of seconds leading up to REM onset (for plotting)
AfterTime = 60;   %desired number of seconds after REM onset
    % for LC use 10 before/after to find peak that ends REM bouts
    %for REM offset definition use 30 before/after to look at baseline levels
    
prism_downsampling = 20;

BeforeSample = round(BeforeTime * signal_fs); %convert time into # of sampling points
AfterSample = round(AfterTime * signal_fs);

% Below I find the time indices of REM onsets - to index into the full time and signal traces
REM_offset_idx = round(REM_offset_cut*signal_fs);

%Extract NE traces before/after REM onset
[REM_offset_LCtraces,ValidLC_event] = TriggerTraces(delta3_filt,REM_offset_idx,BeforeSample,AfterSample);
[REM_offset_NEtraces,ValidNE_event] = TriggerTraces(delta2_filt,REM_offset_idx,BeforeSample,AfterSample);

%downsampled traces
REM_offset_LCtraces_ds = downsample(REM_offset_LCtraces, prism_downsampling);
REM_offset_NEtraces_ds = downsample(REM_offset_NEtraces, prism_downsampling);

%time vector in seconds. Time(0) = REM offset
REMoff_traces_time = ((1:length(REM_offset_LCtraces))/signal_fs) - BeforeTime; 

%plot mean LC and NE traces at REM offset
figure
    plot(REMoff_traces_time, mean(REM_offset_LCtraces, 2), 'DisplayName','LC')
    hold on
    plot(REMoff_traces_time, mean(REM_offset_NEtraces,2), 'DisplayName','NE')
    legend
    xlabel('time (s)');
    ylabel('dF/F');

% plot all NE and LC traces at REM offset
figure
plot(REMoff_traces_time, REM_offset_NEtraces)
figure
plot(REMoff_traces_time, REM_offset_LCtraces)


% mean trace
mean_LC_REMoffset = mean(REM_offset_LCtraces, 2);
mean_LC_REMoffset_ds = downsample(mean_LC_REMoffset, prism_downsampling);
mean_NE_REMoffset = mean(REM_offset_NEtraces,2);
mean_NE_REMoffset_ds = downsample(mean_NE_REMoffset, prism_downsampling);

% SEM of traces
SEM_LC_REMoff = downsample(std(REM_offset_LCtraces, [], 2)./ sqrt(size(REM_offset_LCtraces,2)),prism_downsampling);
SEM_NE_REMoff = downsample(std(REM_offset_NEtraces, [], 2)./ sqrt(size(REM_offset_NEtraces,2)),prism_downsampling);

% input vector for n (for prism)
number_vector_LC_REMoff = zeros(1,length(SEM_LC_REMoff))'+size(REM_offset_LCtraces,2);
number_vector_NE_REMoff = zeros(1,length(SEM_NE_REMoff))'+size(REM_offset_NEtraces,2);

% timeline for prism - starting at 0
prism_timeline = downsample(REMoff_traces_time, prism_downsampling)';
prism_timeline_zero = (prism_timeline-prism_timeline(1));


% REM offset scoring EEG/EMG vs FP scoring
REMoff_diff = zeros(length(REM_offset_cut),1);
for i=1:length(REM_offset_cut)
    REMoff_FPdef = find(REM_offset_NEtraces(:,i)==min(REM_offset_NEtraces(:,i))); %offset defined as the lowest point on NE traces
    %REMoff_FPdef = find(diff(REM_offset_NEtraces(:,i))==max(diff(REM_offset_NEtraces(:,i)))); %offset defined as steapest increase
    REMoff_EEGdef = round(BeforeTime*signal_fs);
    REMoff_diff(i,1) = (REMoff_FPdef-REMoff_EEGdef)/signal_fs; %time (s) of NE minimum relative to EEG offset
end
mean(REMoff_diff);


% Peak characteristics of LC peaks that terminates REM bouts
for i=1:length(REM_offset_cut)
    [LCpks_i, LClocs_i, w_i, p_i] = findpeaks(REM_offset_LCtraces(:,i), REMoff_traces_time, 'MinPeakDistance',2, 'MinPeakProminence', 2);
    figure
        plot(REMoff_traces_time,REM_offset_LCtraces(:,i),'r');
        hold on;
        plot(LClocs_i, LCpks_i,'ko');
        xlabel('time (s)');
    REMoff_LCpks(i,1) = LCpks_i(1); % NB! the first detected peak will be used
    REMoff_LClocs(i,1) = LClocs_i(1);
    REMoff_w(i,1) = w_i(1);
    REMoff_p(i,1) = p_i(1);
    peaktime_idx = find(REMoff_traces_time==LClocs_i(1)); %index of peak time 
    REM_bl_start_idx = peaktime_idx-round(5*signal_fs); % index of timepoint 5 s prior to peak time
    REM_bl_end_idx = peaktime_idx-round(3*signal_fs); % index of timepoint 3 s prior to peak time
    REMoff_basline(i,1) = mean(REM_offset_LCtraces(REM_bl_start_idx:REM_bl_end_idx,i));
    REMoff_LCampl(i,1) = REMoff_LCpks(i,1)-REMoff_basline(i,1);
end


% NE baseline pre REM offset
REM_bl_start_time = 15; %sec prior to transition
REM_bl_end_time = 5; %sec prior to transition
NEbase_preREMoff = zeros(length(REM_offset_cut),1);

for i=1:length(REM_offset_cut)
    REM_bl_start_idx = abs(round(REM_bl_start_time-BeforeTime*signal_fs)); % index of timepoint 5 s prior to peak time
    REM_bl_end_idx = abs(round(REM_bl_end_time-BeforeTime*signal_fs)); % index of timepoint 3 s prior to peak time
    NEbase_preREMoff_i = mean(REM_offset_NEtraces(REM_bl_start_idx:REM_bl_end_idx,i));
    NEbase_preREMoff(i,1) = NEbase_preREMoff_i;
end

% NE peak following REM offset
NEpeak_postREMoff = zeros(length(REM_offset_cut),1);

for i=1:length(REM_offset_cut)
    NEpeak_postREMoff_i = max(REM_offset_NEtraces(BeforeSample:end,i)); %finds the maximum value in the NE trace after REM offset and [AfterTime] seconds after.
    NEpeak_postREMoff(i,1) = NEpeak_postREMoff_i;
end

% NE increase at REM offset
NE_REMoff_increase = NEpeak_postREMoff-NEbase_preREMoff;


% NE rise time at offset
fit_type = fittype('a*exp(-abs(b)*t)+c', 'independent', 't'); % try w/wo c constant and compare decay
f_REMoff = fit(REMoff_traces_time', mean_NE_REMoffset,fit_type); 

time_constant = 1/(f.b); %tau is 1/b - maybe different for rise time
decay_pos = abs(time_constant);

figure
plot(f_REMoff,REMoff_traces_time', mean_NE_REMoffset)


%% 19) LC frequency

minPP = 1.5; 
[LC_pks_full, pklocs_full, w, p] = findpeaks(delta3_filt, sec_signal, 'MinPeakDistance',5, 'MinPeakProminence',minPP);
    % LC_pks = amplitude of LC peaks, pklocs = time of peaks

% plotting the found LC events
figure
    plot(sec_signal(1000:end),delta3_filt(1000:end),'r');
    hold on;
    plot(pklocs_full, delta3_filt(round(pklocs_full*signal_fs)),'ko');
    xlabel('time (s)');
    
% Divide LC event into wake, sws, and REM
[TriggerBool_wake] = IsInInterval(pklocs_full, wake_periods_cut);
[TriggerBool_sws] = IsInInterval(pklocs_full, sws_periods_cut);
[TriggerBool_REM] = IsInInterval(pklocs_full, REM_periods_cut);

%Frequency of LC events
total_wake_duration = sum(wake_duration_cut); % total duration in seconds
LCfreq_wake = sum(TriggerBool_wake)/total_wake_duration; % Hz
LCfreq_wake = sum(TriggerBool_wake)/(total_wake_duration/60); % per minute

total_sws_duration = sum(sws_duration_cut); % total duration in seconds
LCfreq_sws = sum(TriggerBool_sws)/total_sws_duration; % Hz
LCfreq_sws = sum(TriggerBool_sws)/(total_sws_duration/60); % per minute

total_REM_duration = sum(REM_duration_cut); % total duration in seconds
LCfreq_REM = sum(TriggerBool_REM)/total_REM_duration; % Hz
LCfreq_REM = sum(TriggerBool_REM)/(total_REM_duration/60); % per minute


%% 20) Manual REM score

REMon_manual = EEG_FP_triple_237_manualREMscore(:,1);
REMoff_manual = EEG_FP_triple_237_manualREMscore(:,2);

BeforeTime = 60;  %desired number of seconds leading up to
AfterTime = 60;   %desired number of seconds after

BeforeSample = round(BeforeTime * sampling_freq); %convert time into # of sampling points
AfterSample = round(AfterTime * sampling_freq);

% Below I find the time indices of REM onsets - to index into
REMon_manual_idx = round(REMon_manual*sampling_freq);
REMoff_manual_idx = round(REMoff_manual*sampling_freq);
%REMon_EEG_idx = round(REM_onset_cut*sampling_freq); %EEG score
% REMoff_EEG_idx = round(REM_offset_cut*sampling_freq); %EEG score

%Extract NE and LC traces before/after REM onset
[REM_onset_EEGtraces,ValidLC_event] = TriggerTraces(EEG_rawtrace_cut,REMon_EEG_idx,BeforeSample,AfterSample);
[REM_offset_EEGtraces,ValidNE_event] = TriggerTraces(EEG_rawtrace_cut,REMoff_EEG_idx,BeforeSample,AfterSample);

REMon_EEG_time = ((1:length(REM_onset_EEGtraces))/sampling_freq) - BeforeTime; %time vector in seconds. Time(0) = LC event 

%plot mean LC and NE traces at REM onset
figure
    a = subplot(2, 1, 1);
    plot(REMon_EEG_time, mean(REM_onset_EEGtraces, 2), 'DisplayName','REM onset')
    hold on
    b = subplot (2, 1, 2);
    plot(REMon_EEG_time, mean(REM_offset_EEGtraces,2), 'DisplayName','REM offset')
    legend
    xlabel('time (s)');
    ylabel('dF/F');



% plot all NE traces triggered by LC
figure
plot(REMon_EEG_time, REM_onset_EEGtraces)
figure
plot(REMon_EEG_time, REM_offset_EEGtraces)


%% 21) EEG power spectrum analysis
fs_signal = 1:1:length(delta2_filt);
sec_signal = fs_signal/signal_fs;

window = 5; %sec. 1 for 30 sec

power_bands = {[1, 4], [4, 8], [8, 15], [15, 30]}; % delta was 0.2 before
%power_bands = {[0.2, 4], [4, 7], [7, 15], [15, 30] [7 10]};
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
    %plot( downsample(sec_signal, 100)+mouse{6},downsample(delta2_filt, 100))
    plot( downsample(sec_signal, 100),downsample(delta2_filt, 100))
 
    b = subplot(4, 1, 2);
    imagesc(time_spectrogram_zero, F, filtered_mean_spectrogram); %plot the log spectrum
    set(gca,'YDir', 'normal'); % flip the Y Axis so lower frequencies are at the bottom
    ylim([0, 30]);
    caxis([-6.7, -4])
    %     colorbar()
    %h = colorbar;
    %title('#237');
    colormap(gca, 'parula');
    hold on
    for band_i = 1:length(power_bands)
        plot([-295, -295], power_bands{band_i}, 'LineWidth', 5)
    end
    
    c = subplot(4, 1, 3);
    band_power_collector = [T];
    
   %[val,idx]=min(abs(T-norm_time));
    %minVal =T(idx);
    
    %[val1, idx1] = min(abs(F-power_band(1)));
    %[val2, idx2] = min(abs(F-power_band(2)));
    %normalization_factor = mean(mean(mean_spectrogram((idx1:idx2)), 1:find(T==200)));
    %normalization_factor = mean(mean(mean_spectrogram(find(F==total_power_band(1)):find(F==total_power_band(2)), 1:find(T==minVal))));
    normalization_factor = mean(mean(mean_spectrogram(find(F==total_power_band(1)):find(F==total_power_band(2)), 1:find(T==round(1)))));
    
    for band_i = 1:length(power_bands)
        power_band = power_bands{band_i};
       %[val1, idx1] = min(abs(F-power_band(1)));
       %[val2, idx2] = min(abs(F-power_band(2)));
       % power_trace = mean(mean_spectrogram((idx1:idx2)), 1);
        power_trace = mean(mean_spectrogram(find(F==power_band(1)):find(F==power_band(2)), :), 1);
        %normalized_power_trace = power_trace/-normalization_factor+2;
        normalized_power_trace = power_trace;
        band_power_collector = [band_power_collector; normalized_power_trace];
        plot(time_spectrogram_zero, normalized_power_trace)
        hold on
    end
    d = subplot(4, 1, 4);
    sigma = band_power_collector(4,:); %theta/delra ratio
    plot(time_spectrogram_zero, sigma)
    title ('sigma')
  
    linkaxes([a,b,c,d],'x');
    
%% 22) Plot NE, EEG powerspectrum, and sigma
% And select time points for sigma extraction

% NE manual scores
EEG_FP_237_NE_manualscore = 'Q:\Personal_folders\Celia Kjaerby\Fiber photometry\20200120_EEG_FP\20200116_237_EEG_FP\237 NE peak-trough times (manual score).xlsx';
EEG_FP_243_NE_manualscore = 'Q:\Personal_folders\Celia Kjaerby\Fiber photometry\20200120_EEG_FP\20200123_243_EEG_FP\243 NE peak-trough times (manual score).xlsx';
EEG_FP_319_NE_manualscore = 'Q:\Personal_folders\Celia Kjaerby\Fiber photometry\20200629_EEG_FP_batchII\20200622_319_FP_EEG\20200622_319_LC_NE_As\319 NE peak-trough times (manual score).xlsx';
EEG_FP_331_NE_manualscore = 'Q:\Personal_folders\Celia Kjaerby\Fiber photometry\20200629_EEG_FP_batchII\20200624_331_FP_EEG\20200624_331_LC_NE_As\331 NE peak-trough times (manual score).xlsx';
EEG_FP_307_NE_manualscore = 'Q:\Personal_folders\Celia Kjaerby\Fiber photometry\20200629_EEG_FP_batchII\20200626_307_FP_EEG\20200626_307_LC_NE_As\307 NE peak-trough times (manual score).xlsx';
EEG_FP_323_NE_manualscore = 'Q:\Personal_folders\Celia Kjaerby\Fiber photometry\20200629_EEG_FP_batchII\20200625_323_FP_EEG\20200625_323_LC_NE_As\323 NE peak-trough times (manual score).xlsx';
EEG_FP_329_NE_manualscore = 'Q:\Personal_folders\Celia Kjaerby\Fiber photometry\20200629_EEG_FP_batchII\20200624_329_FP_EEG\20200624_329_LC_NE_As\329 NE peak-trough times (manual score).xlsx';

manualscore = EEG_FP_329_NE_manualscore;

manual_NREM_peak = xlsread(manualscore, 'NREM', 'A:A');
manual_MA_peak = xlsread(manualscore, 'MA', 'A:A');
manual_wake_peak = xlsread(manualscore, 'wake', 'A:A');
manual_REM_peak = xlsread(manualscore, 'REM', 'A:A');

Time_points = manual_NREM_peak;

figure
a = subplot(4,1,1);
    plot(ds_sec_signal, ds_delta2_filt)
    hold on
        plot(Time_points, delta2_filt(round(Time_points*signal_fs)), 'r*')
    xlabel('time (s)');
    ylabel('dF/F (%)');
    title('NE2m');
b = subplot(4,1,2);
    imagesc(time_spectrogram_zero, F, filtered_mean_spectrogram); %plot the log spectrum
    set(gca,'YDir', 'normal'); % flip the Y Axis so lower frequencies are at the bottom
    ylim([0, 30]);
    caxis([-6.7, -4])
    %     colorbar()
    %h = colorbar;
    %title('#237');
    colormap(gca, 'parula');
    hold on
    for band_i = 1:length(power_bands)
        plot([-295, -295], power_bands{band_i}, 'LineWidth', 5)
    end
c = subplot(4,1,3);
    plot(time_spectrogram_zero, sigma)
    xlabel('time (s)');
    ylabel('power');
    title('sigma');
d = subplot(4,1,4);
    ds_EEG_rawtrace = downsample(EEG_rawtrace_cut, 10);
    plot_sleep(ds_EEG_time, ds_EEG_rawtrace, sleepscore_time_cut, wake_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut);
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([a,b,c,d,e],'x');

%% 23) EEG sigma extraction

Time_points = EEG_FP_triple_237_NREM_NE_onset;

time_before = 65; %sec
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
    %if (sound_time+time_after)*signal_fs < length(delta2_filt)
        sound_index = round(sound_time*frq);
        %sound_index_FP = round(sound_time*signal_fs);
        sound_before_index = round(sound_index-frq*time_before);
        sound_after_index = round(sound_index+frq*time_after);
        %sound_before_index_FP = round(sound_index_FP-signal_fs*time_before);
        %sound_after_index_FP = round(sound_index_FP+signal_fs*time_after);
        %fp_NE_trace = delta2_filt(sound_before_index_FP:sound_after_index_FP);
        eeg_sound_trace = EEG_rawtrace_cut(:, sound_before_index:sound_after_index);
        [sound_spectrogram, F, T] = spectrogram(eeg_sound_trace,round(frq*analysis_window),[],frw,frq,'yaxis'); % F = frequenciy vector, T=time vector
        sigma_power_trace = mean(sound_spectrogram(find(F==8):find(F==15), :), 1);
        normalized_sigma_power = sigma_power_trace/-normalization_factor+2;
        sigma_epocs = [sigma_epocs  normalized_sigma_power'];
        sound_spectrogram_collector = cat(3, sound_spectrogram_collector, sound_spectrogram);
        %sound_FP_NE_collector = [sound_FP_NE_collector fp_NE_trace'];
    %else continue
    %end
end

mean_spectrogram = nanmean(log(abs(sound_spectrogram_collector)), 3);

time_spectrogram_zero = T-time_before; % to get REM onset at 0 instead of 300
%time_FP = (1:1:length(fp_NE_trace))/signal_fs -time_before;

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
 
sigma_max_20pre_sound = prctile(sigma_mean(round((time_before-20)*specto_fs:time_before*specto_fs)),95); % top 95 percentile of sigma trace 20 s pre sound
sigma_min_20post_sound = prctile(sigma_mean(round(time_before*specto_fs:(time_before+20)*specto_fs)),5); % bottom 5 percentile of sigma trace 20 s post sound


%% 24) Re-classify MAs as NREM sleep (for NE oscillation rate)

% Creating one vector with different behaviors represented by unique
% numbers (1=wake, 4=sws, 9=REM, 15=MA) at frequency 1Hz
boutscore_vector = zeros([1, (Info.HypnoFiles.Duration)+time_correction+6]);

% Here using the unaligned "uncut" vectors
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

for i=1:length(REM_onset)
    t = REM_onset(i)+1;
    d = REM_duration(i)-1;
    boutscore_vector(t:t+d) = 9; %REM=9
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
%NREMinclMA_duration = NREMinclMA_duration(NREMinclMA_duration>14); % excluding MAs during light/desipramine sleep
%NREMinclMA_onset = NREMinclMA_onset(NREMinclMA_duration>14);
%NREMinclMA_offset = NREMinclMA_offset(NREMinclMA_duration>14);

%% 25) NE dynamic range

wake_woMA_periods_cut_rmREMtail = wake_woMA_periods_cut;
wake_woMA_periods_cut_rmREMtail(1,1)=1; % remove first second or trace
for i = 1:length(wake_woMA_periods_cut)
    if any(wake_woMA_periods_cut(i,1)==REM_periods_cut(:,2)')   % remove 20 s of wake bouts follwoing REM bouts to exclude low transition values
        wake_woMA_periods_cut_rmREMtail(i,1)= wake_woMA_periods_cut(i,1)+20;
    else
        continue
    end
end

wake_woMA_periods_cut_rmREMtail_durations = wake_woMA_periods_cut_rmREMtail(:,2)-wake_woMA_periods_cut_rmREMtail(:,1);
wake_woMA_periods_cut_rmREMtail_filt = wake_woMA_periods_cut_rmREMtail(wake_woMA_periods_cut_rmREMtail_durations>0,:);

delta2_filt_positive = delta2_filt+abs(min(delta2_filt(1000:end)));

% divide NE trace based on brain states
REM_periods_cut_sampl = REM_periods_cut*signal_fs; %converts on-/offset times to sampling points in FP dataset
NREMinclMA_periods_cut_sampl = NREMinclMA_periods_cut*signal_fs;
wake_woMA_periods_cut_sampl = wake_woMA_periods_cut_rmREMtail_filt*signal_fs;

Delta2_REM = [];
NE_REM_CoVa = [];
NE_REM_rms = [];
for i = 1:size(REM_periods_cut_sampl, 1)
   period1 = REM_periods_cut_sampl(i,:);
   if period1(2)<length(delta2_filt_positive)
   trace2 = delta2_filt_positive(period1(1):period1(2));
   CoVa =(std(trace2)/mean(trace2))*100;
   rootms = rms(trace2);
   NE_REM_CoVa = [NE_REM_CoVa CoVa];
   NE_REM_rms = [NE_REM_rms rootms];
   Delta2_REM = [Delta2_REM trace2];
   else
       continue
   end
end

Delta2_sws = [];
NE_sws_CoVa = [];
NE_sws_rms = [];
for i = 1:size(NREMinclMA_periods_cut_sampl,1)
   period1 = NREMinclMA_periods_cut_sampl(i,:)+1;
   if period1(2)<length(delta2_filt_positive)
   trace2 = delta2_filt_positive(period1(1):period1(2));
   CoVa =(std(trace2)/mean(trace2))*100;
   rootms = rms(trace2);
   NE_sws_CoVa = [NE_sws_CoVa CoVa];
   NE_sws_rms = [NE_sws_rms rootms];
   Delta2_sws = [Delta2_sws trace2];
   else
       continue
   end
end

Delta2_wake = [];
NE_wake_CoVa = [];
NE_wake_rms = [];
for i = 1:size(wake_woMA_periods_cut_sampl,1)
   period1 = wake_woMA_periods_cut_sampl(i,:);
   if period1(2)<length(delta2_filt_positive)
   trace2 = delta2_filt_positive(period1(1)+1:period1(2));
   CoVa =(std(trace2)/mean(trace2))*100;
   rootms = rms(trace2);
   NE_wake_CoVa = [NE_wake_CoVa CoVa];
   NE_wake_rms = [NE_wake_rms rootms];
   Delta2_wake = [Delta2_wake trace2];
   else
       continue
   end
end

figure
a = subplot(3,1,1);
plot(Delta2_sws);
title 'sws'
b = subplot(3,1,2);
plot(Delta2_wake);
title 'wake'
c = subplot(3,1,3);
plot(Delta2_REM);
title 'REM'

%deffine the edges to cover the values of the histrogram data. 
edges = [-20:0.5:20];

% histogram plots
figure
a = subplot(3,1,1);
Histo_sws = histogram(Delta2_sws, edges);
title 'sws'
b = subplot(3,1,2);
Histo_wake = histogram(Delta2_wake, edges);
title 'wake'
c = subplot(3,1,3);
Histo_REM = histogram(Delta2_REM, edges);
title 'REM'

Histo = [];
histo_time = Histo_REM.BinEdges';
Histo(:,1) = histo_time(2:end);
Histo(:,2) = Histo_wake.BinCounts';
Histo(:,3) = Histo_sws.BinCounts';
Histo(:,4) = Histo_REM.BinCounts';

NE_sws_CoVa_all =(std(Delta2_sws)/mean(Delta2_sws))*100;
NE_wake_CoVa_all =(std(Delta2_wake)/mean(Delta2_wake))*100;
NE_REM_CoVa_all =(std(Delta2_REM)/mean(Delta2_REM))*100;

NE_sws_CoVa_mean = mean(NE_sws_CoVa);
NE_wake_CoVa_mean = mean(NE_wake_CoVa);
NE_REM_CoVa_mean = mean(NE_REM_CoVa);

NE_REM_rms_mean = mean(NE_REM_rms);
NE_sws_rms_mean = mean(NE_sws_rms);
NE_wake_rms_mean = mean(NE_wake_rms);

NE_sws_SEM =(std(Delta2_sws)/sqrt(length(Delta2_sws)));
NE_wake_SEM =(std(Delta2_wake)/sqrt(length(Delta2_wake)));
NE_REM_SEM =(std(Delta2_REM)/sqrt(length(Delta2_REM)));

NE_sws_var = var(Delta2_sws);
NE_wake_var = var(Delta2_wake);
NE_REM_var = var(Delta2_REM);

%% NE oscillation PSD

psd_periods = wake_woMA_periods_cut_sampl;
data_bouts = cell(1, length(psd_periods));

fc = 0.025; % low cutoff filter
filter_fs = 10; 

PXX_psd = [];

for i=1:length(psd_periods)%-1
    data_bouts{i} = delta2_filt(psd_periods(i,1):psd_periods(i,2));
    ds_data_bout_i = downsample(data_bouts{i},100);
    [b,a] = butter(100,fc); % 100th order filter
    dataOut = filter(b,a,ds_data_bout_i);
    [pxx_psd, f_psd] = pwelch(dataOut, [], [],[0:0.002:1], 10);
    logpxx_psd = 10*log10(pxx_psd);
    FX_psd{i} = f_psd;
    PXX_psd(:,i) = logpxx_psd;
end

mean_PXX_psd = mean(PXX_psd,2);

prism_psd_psd = mean_PXX_psd(f_psd<1);
prism_freq = f_psd(f_psd<1);

% power spectral density plot
figure
plot(prism_freq, prism_psd_psd)




[b,a] = butter(100,fc/(filter_fs/2));
dataIn = delta2_filt;
dataOut = filter(b,a,dataIn);

%% NE FFT (fast fourier transform)

psd_periods = NREMinclMA_periods_cut_sampl;
data_bouts = cell(1, length(psd_periods));

PXX_psd = [];
psd_data_collect = [];

for i=1:3%length(psd_periods)%-1
    data_bouts{i} = delta2_filt(psd_periods(i,1):psd_periods(i,2));
    fft_data = fft(data_bouts{i});
    P2 = abs(fft_data/length(data_bouts{i}));
    P1 = P2(1:length(data_bouts{i})/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = signal_fs*(0:(length(data_bouts{i})/2))/length(data_bouts{i});
    figure
    plot(f,P1) 
    title('Single-Sided Amplitude Spectrum of X(t)')
    xlim([0 1])
    xlabel('f (Hz)')
    ylabel('|P1(f)|')
end










