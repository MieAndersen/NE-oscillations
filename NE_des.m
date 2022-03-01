%% 1) Define mouse data

% data structure:
    % 1) FP raw data
    % 2) EEG raw data
    % 3) EEG sleep score
    % 4) 465 channel name
    % 5) 405 channel name (only batch II)
    % 6) time of injection (s) - according to EEG/EMG recording
    % 7) Hypnogram time correction (s)
    % 8) Interval used for signal normalization (only batch II)

example_mouseID = {'C:\Users\username\data\FP_data_folder' 'C:\Users\username\data\EEG_data.exp' 'C:\Users\username\data\sleep_score.xlsx' 'channel 465' 'channel 405' 2000 0 (1000:1800)};

mouse = example_mouseID;

%% 2) Load FP data

data = TDTbin2mat(mouse{1});  % custom function provided by Tucker Davis Technologies

signal_fs = data.streams.(mouse{4}).fs;
signal_465 = data.streams.(mouse{4}).data; %hSyn-NE

% removing FP trace prior to first TTL pulse
TTL_FP = data.epocs.PtC0.onset;
TTL_gap = diff(TTL_FP) > 5 + 1;
if isempty(find(TTL_gap == 1, 1))
    TTL_onset = TTL_FP(1);
else 
    TTL_onset = TTL_FP(find(TTL_gap==1)+1);
end

first_TTL = TTL_onset(1)*signal_fs;
onset_FP = first_TTL;

signal_465 = signal_465(round(onset_FP):end);

%% 3a) Normalize and plot (batch I)
% Normalize using median

fs_signal = 1:length(signal_465);
sec_signal = fs_signal/signal_fs;

% deltaF/F
med_465 = median(signal_465(1:round(30*60*signal_fs))); %median of 30 minute baseline is used
delta_465 = ((signal_465 - med_465)/med_465)*100;
delta465_filt = filtfilt(MeanFilter,1,double(delta_465));
ds_delta465_filt = downsample(delta465_filt, 100);

fs_signal = 1:1:length(delta465_filt);
sec_signal = fs_signal/signal_fs;
ds_sec_signal = downsample(sec_signal, 100);

figure
plot(ds_sec_signal(1000:end), ds_delta465_filt(1000:end))
title('NE2m');


%% 3b) Normalize and plot (batch II)
% Normalize using 405 nm channel

signal_405 = data.streams.(mouse{5}).data; %autofluorescence
signal_405 = signal_405(onset_FP:end);

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
ds_sec_signal = downsample(sec_signal, ds_factor_FP);

figure
plot(ds_sec_signal, ds_delta465_filt)
title('NE2m');

%% 4) loading and plotting EEG and EMG raw data

Info=loadEXP([mouse{2}],'no'); % custom function provided by Viewpoint Behavior Technology

TimeReldebSec=0; %start extract data from the beginning (first bin)
TimeRelEndSec=Info.BinFiles.Duration; %inf to include all data (including last bin)

[Data,Time]=ExtractContinuousData([],Info,[],TimeReldebSec, TimeRelEndSec,[],1);

EMG_rawtrace = Data(1,1:end);
EEG_rawtrace = Data(2,1:end);

sampling_freq = Info.Fs;
EEG_time = (1:length(EEG_rawtrace))/sampling_freq;

% Plot of EEG and EMG traces
figure
h(1) = subplot(2,1,1);
    plot(EEG_time, EMG_rawtrace); 
    xlabel('time (s)');
    ylabel('EMG (V)');
h(2) = subplot(2,1,2);
    plot(EEG_time, EEG_rawtrace); 
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([h(1),h(2)],'x');

%% 5) open EEG scoring

time_correction = mouse{7}; % NB! If there is a systematic time lag between EEG/EMG traces and scoring adjust for it by seconds here
EEG_sleepscore = xlsread(mouse{3});

%Awake
wake_onset = rmmissing(EEG_sleepscore(:, 2));
wake_duration = rmmissing(EEG_sleepscore(:, 3)); 

%Slow-wave sleep
sws_onset = rmmissing(EEG_sleepscore(:, 6)); 
duration_sws = rmmissing(EEG_sleepscore(:, 7)); 

%undefined sleep due to des
des_onset = rmmissing(EEG_sleepscore(:, 10)); 
des_duration = rmmissing(EEG_sleepscore(:, 11)); 

%REM
if size(EEG_sleepscore,2) < 13 % in case of no REM bouts
    REM_onset = NaN;
    REM_duration = NaN;
else
    REM_onset = rmmissing(EEG_sleepscore(:, 14)); 
    REM_duration = rmmissing(EEG_sleepscore(:, 15)); 
end

% Most EEG scorings don't start at time 0 which shifts the timeline of the
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

% Create binary vectors for sleep stages
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

des_binary_vector = zeros([1, (Info.HypnoFiles.Duration)+6]); 
for i=1:length(des_onset)
    t = des_onset(i)+1;
    d = des_duration(i)-1;
    des_binary_vector(t:t+d) = 1;
end

sleepscore_time = 0:length(wake_binary_vector)-1; % Should be same length for wake/sws/REM binary vectors

% check alignment of scoring and adjust time_correction if necessary
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

wake_periods = [wake_onset wake_onset+wake_duration];
sws_periods = [sws_onset sws_onset+duration_sws];
REM_periods = [REM_onset REM_onset+REM_duration];
des_periods = [des_onset des_onset+des_duration];

%% 6) Dividing wake bouts into microarousals (MA) and wake w/o MA

MA_maxdur = 15; % maximum duration of microarrousal
MA_idx = find(wake_duration < MA_maxdur);
MA_onset = wake_onset(MA_idx);
MA_duration = wake_duration(MA_idx);
MA_binary_vector = zeros([1, (Info.HypnoFiles.Duration)+6]); 
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
wake_woMA_binary_vector = zeros([1, (Info.HypnoFiles.Duration)+6]); 
for i=1:length(wake_woMA_onset) 
    t = wake_woMA_onset(i)+1;
    d = wake_woMA_duration(i)-1;
    wake_woMA_binary_vector(t:t+d) = 1;
end

% 2-column vectors with on- and offsets for each state
MA_periods = [MA_onset MA_onset+MA_duration];
wake_woMA_periods = [wake_woMA_onset wake_woMA_onset+wake_woMA_duration];


%% 7) Alingment of EEG recording and FP recording

% TTL pulse from FP
TTL_pulse = Data(3,1:end);
onset_EEG = find(diff(TTL_pulse>1*10^-3));
onset_EEG_time = onset_EEG/sampling_freq;
onset_EEG_time_diff = diff(onset_EEG_time);

TTL_gap_EEG = onset_EEG_time_diff > 5;
if isempty(find(TTL_gap_EEG==1, 1))
    onset_EEG = onset_EEG(1);
else 
    onset_EEG = onset_EEG(find(onset_EEG_time_diff>5)+1);
end

TTL_EEG_onset = onset_EEG/sampling_freq;

%Cutting EEG/EMG traces leading up to first TTL 
EMG_rawtrace_cut = EMG_rawtrace(round(TTL_EEG_onset*sampling_freq):end);
EEG_rawtrace_cut = EEG_rawtrace(round(TTL_EEG_onset*sampling_freq):end);
EEG_time_cut = (1:length(EEG_rawtrace_cut))/sampling_freq;

% Remove time before TTL from EEG score to align with FP trace
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

if ~isnan(REM_onset)
    [REM_onset_cut, REM_offset_cut] = binary_to_OnOff(REM_binary_vector_cut);
    REM_duration_cut = REM_offset_cut - REM_onset_cut;
else
    REM_onset_cut = NaN;
    REM_offset_cut = NaN;
    REM_duration_cut = NaN;
end

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

% Alignement of transitions vector:
transition_sws_wake_cut = transition_sws_wake-TTL_EEG_onset;
transition_sws_MA_cut = transition_sws_MA-TTL_EEG_onset;


%% 8) Plotting all traces and scorings together 

sleepscore_time_cut = 0:length(wake_binary_vector_cut)-1; % should be same length for wake/sws/REM

fig = figure; 
a = subplot(3,1,1);
    plot_sleep(ds_sec_signal, ds_delta465_filt, sleepscore_time_cut, wake_woMA_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut, MA_binary_vector_cut);
    title('NE2m');
b = subplot (3,1,2);
    ds_EEG_time = downsample(EEG_time_cut, 10);
    ds_EMG_rawtrace = downsample(EMG_rawtrace_cut, 10);
    plot_sleep(ds_EEG_time, ds_EMG_rawtrace, sleepscore_time_cut, wake_woMA_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut, MA_binary_vector_cut);
    xlabel('time (s)');
    ylabel('EMG (V)');
c = subplot(3,1,3);
    ds_EEG_rawtrace = downsample(EEG_rawtrace_cut, 10);
    plot_sleep(ds_EEG_time, ds_EEG_rawtrace, sleepscore_time_cut, wake_woMA_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut, MA_binary_vector_cut);
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([a,b,c],'x');


%% 9) Sleep quantification and PSD

injection_time = mouse{6};
analysis_start_time = mouse{6}+ 60*60; % analysis start 1 h after injection
analysis_end_time = analysis_start_time + 2.5*60*60; % analyze 2.5 hours of recording

sws_onset_exclBL = sws_onset(sws_onset>analysis_start_time & sws_periods(:,2)<(analysis_end_time));
duration_sws_exclBL = duration_sws(sws_onset>analysis_start_time & sws_periods(:,2)<(analysis_end_time));

REM_onset_exclBL = REM_onset(REM_onset>analysis_start_time & REM_periods(:,2)<(analysis_end_time));
REM_duration_exclBL = REM_duration(REM_onset>analysis_start_time & REM_periods(:,2)<(analysis_end_time));

MA_onset_exclBL = MA_onset(MA_onset>analysis_start_time & MA_periods(:,2)<(analysis_end_time));
MA_duration_exclBL = MA_duration(MA_onset>analysis_start_time & MA_periods(:,2)<(analysis_end_time));

wake_woMA_onset_exclBL = wake_woMA_onset(wake_woMA_onset>analysis_start_time & wake_woMA_periods(:,2)<(analysis_end_time));
wake_woMA_duration_exclBL = wake_woMA_duration(wake_woMA_onset>analysis_start_time & wake_woMA_periods(:,2)<(analysis_end_time));

des_onset_exclBL = des_onset(des_onset>analysis_start_time & des_periods(:,2)<(analysis_end_time));
des_duration_exclBL = des_duration(des_onset>analysis_start_time & des_periods(:,2)<(analysis_end_time));

prism_exp_totaltime = [sum(wake_woMA_duration_exclBL)+ 2470, sum(duration_sws_exclBL), sum(REM_duration_exclBL), sum(MA_duration_exclBL), sum(des_duration_exclBL)];
excl_flank_dur = (2*60*60)-sum(prism_exp_totaltime); 
prism_exp_minperhour = prism_exp_totaltime/(sum(prism_exp_totaltime)/60);
prism_exp_secperhour = prism_exp_minperhour*60;
prism_exp_n = [length(wake_woMA_duration_exclBL), length(duration_sws_exclBL), length(REM_duration_exclBL), length(MA_duration_exclBL)];
prism_exp_n_perhour = prism_exp_n/(sum(prism_exp_totaltime)/60/60); 
prism_exp_boutduration = [mean(wake_woMA_duration_exclBL), mean(duration_sws_exclBL), mean(REM_duration_exclBL), mean(MA_duration_exclBL)];
prism_MA_per_NREMdur = prism_exp_n(4)/(prism_exp_totaltime(2)/60/60); 


% Power spectral densities
t1 = sws_onset_exclBL;
t2 = sws_onset_exclBL+duration_sws_exclBL;

tsamp1 = floor(t1*sampling_freq); %eeg start time 
tsamp2 = floor(t2*sampling_freq); %eeg end time
NREM_data = cell(1, numel(tsamp1));

PXX = [];
NREM_data_collect = [];

for i=1:numel(tsamp1)
    NREM_data{i} = EEG_rawtrace(tsamp1(i):tsamp2(i));
    NREM_data_cut = EEG_rawtrace(tsamp1(i):tsamp2(i));
    NREM_data_collect = [NREM_data_collect NREM_data_cut];
    [pxx, f] = pwelch(NREM_data{i}, [], [],[0:0.2:100], sampling_freq);
    logpxx = 10*log10(pxx);
    FX{i} = f;
    PXX(:,i) = logpxx;
    PXX(:,i) = pxx;
end

mean_PXX = mean(PXX,2);

prism_psd = mean_PXX(f<45);
prism_freq = f(f<45);

% PSD plot
figure
plot(prism_freq,prism_psd)

prism_mean_sigma_power_density = mean(mean_PXX(f>8 & f<15));


%% 10) mean NE during sleep wake

% select bouts for analysis
indices_sws = find(sws_periods_cut(:,1)>analysis_start_time & sws_periods_cut(:,1)< analysis_end_time);
sws_periods_cut1 = sws_periods_cut(indices_sws,:);

indices_wake = find(wake_periods_cut(:,1)>analysis_start_time & wake_periods_cut(:,1)< analysis_end_time);
wake_periods_cut1 = wake_periods_cut(indices_wake,:);

indices_rem = find(REM_periods_cut(:,1)>analysis_start_time & REM_periods_cut(:,1)< analysis_end_time);
if isempty(indices_rem)
   REM_periods_cut1 = [];
else
   REM_periods_cut1 = REM_periods_cut(indices_rem,:);
end

wake_onset_cut1 = wake_periods_cut1(:,1);
wake_offset_cut1 = wake_periods_cut1(:,2);

sws_onset_cut1 = sws_periods_cut1(:,1);
sws_offset_cut1 = sws_periods_cut1(:,2);

if isempty(REM_periods_cut1)
    REM_onset_cut1 = [];
    REM_offset_cut1 = [];
else
    REM_onset_cut1 = REM_periods_cut1(:,1);
    REM_offset_cut1 = REM_periods_cut1(:,2);
end

REM_periods_cut_fs = REM_periods_cut1*signal_fs;
sws_periods_cut_fs = sws_periods_cut1*signal_fs;
wake_periods_cut_fs = wake_periods_cut1*signal_fs;

% mean NE
Delta2_REM = [];
for i = 1:length(REM_periods_cut_fs)
   period1 = REM_periods_cut_fs(i,:);
   trace2 = delta465_filt(period1(1):period1(2));
   Delta2_REM = [Delta2_REM trace2];
end

Delta2_sws = [];
for i = 1:length(sws_periods_cut_fs)
   period1 = sws_periods_cut_fs(i,:);
   trace2 = delta465_filt(period1(1):period1(2));
   Delta2_sws = [Delta2_sws trace2];
end

Delta2_wake = [];
for i = 1:length(wake_periods_cut_fs)
   period1 = wake_periods_cut_fs(i,:);
   trace2 = delta465_filt(period1(1)+1:period1(2));
   Delta2_wake = [Delta2_wake trace2];
end


ds_delta_sws = downsample(Delta2_sws, 1000)';
ds_delta_wake = downsample(Delta2_wake, 1000)';
ds_delta_REM = downsample(Delta2_REM, 1000)';

mean_sws = mean(ds_delta_sws);
mean_wake = mean(ds_delta_wake);
mean_REM = mean(ds_delta_REM);

prism_mean_NE = [mean_wake  mean_sws mean_REM]';


%% 11) EEG power spectrum analysis

analysis_window = 5;
power_bands = {[1, 4], [4, 8], [8, 15], [15, 30]};
total_power_band = [0, 30];
frw = 0:0.2:30;

Data_EEG = EEG_rawtrace_cut;
frq = sampling_freq;

[transition_spectrogram, F, T] = spectrogram(Data_EEG,round(frq*analysis_window),[],frw,frq,'yaxis');
mean_spectrogram = log(abs(transition_spectrogram));
time_spectrogram_zero = T; 
filtered_mean_spectrogram = imgaussfilt(mean_spectrogram, 4);

specto_fs = length(T)/T(end);

band_power_collector = [T];

normalization_factor = mean(mean(mean_spectrogram(find(F==total_power_band(1)):find(F==total_power_band(2)), 1:find(T==round(1)))));

for band_i = 1:length(power_bands)
    power_band = power_bands{band_i};
    power_trace = mean(filtered_mean_spectrogram(find(F==power_band(1)):find(F==power_band(2)), :), 1);
    normalized_power_trace = power_trace;
    band_power_collector = [band_power_collector; normalized_power_trace];
    plot(time_spectrogram_zero, normalized_power_trace)
    hold on
end

sigma = band_power_collector(4,:);
delta = band_power_collector(2,:);
theta = band_power_collector(3,:);
beta = band_power_collector(5,:);


%% 12) EEG sigma epoc extraction

% transition from NREM to wake or MA within analysis time frame
transition_sws_wake_analysis = transition_sws_wake_cut(transition_sws_wake_cut > injection_time & transition_sws_wake_cut < analysis_end_time); 
transition_sws_MA_analysis = transition_sws_MA_cut(transition_sws_MA_cut > injection_time & transition_sws_MA_cut < analysis_end_time); 

Time_points = transition_sws_wake_analysis;

time_before = 60; %sec
time_after = 60; %sec 
 
analysis_window = 5; %sec. 1 for 30 sec

power_bands = {[1, 4], [4, 8], [8, 15], [15, 30] [6 8]}; 
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
        [sound_spectrogram, F, T] = spectrogram(eeg_sound_trace,round(frq*analysis_window),[],frw,frq,'yaxis');
        sigma_power_trace = mean(sound_spectrogram(find(F==8):find(F==15), :), 1);
        normalized_sigma_power = sigma_power_trace/-normalization_factor+2;
        sigma_epocs = [sigma_epocs  normalized_sigma_power'];
        sound_spectrogram_collector = cat(3, sound_spectrogram_collector, sound_spectrogram);
        sound_FP_NE_collector = [sound_FP_NE_collector fp_NE_trace'];
    else continue
    end
end

mean_spectrogram = nanmean(log(abs(sound_spectrogram_collector)), 3);
time_spectrogram_zero = T-time_before; 
time_FP = (1:1:length(fp_NE_trace))/signal_fs -time_before;

figure()
a = subplot(4, 1, 1);
    filtered_mean_spectrogram = imgaussfilt(mean_spectrogram, 4);
    imagesc(time_spectrogram_zero, F, filtered_mean_spectrogram); 
    set(gca,'YDir', 'normal'); 
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
 
theta_mean = band_power_collector(3,:);

sigma_max_20pre_sound = prctile(sigma_mean(round((time_before-20)*specto_fs:time_before*specto_fs)),95); 
sigma_min_20post_sound = prctile(sigma_mean(round(time_before*specto_fs:(time_before+20)*specto_fs)),5); 
sigma_diff = sigma_max_20pre_sound-sigma_min_20post_sound;
sigma_sum = [sigma_max_20pre_sound sigma_min_20post_sound sigma_diff]';

theta_max_20pre_sound = prctile(theta_mean(round((time_before-20)*specto_fs:time_before*specto_fs)),95); 
theta_min_20post_sound = prctile(theta_mean(round(time_before*specto_fs:(time_before+20)*specto_fs)),5); 
theta_diff = theta_max_20pre_sound-theta_min_20post_sound;
theta_sum = [theta_max_20pre_sound theta_min_20post_sound theta_diff]';

figure
  plot(time_spectrogram_zero, theta_mean)

  
%% 13) NE amplitude and oscillation rate

period = wake_periods_cut1;
period_fs = period*signal_fs;

stitch = [];
trace_2_stitch = [];
for i = 1:length(period_fs)
   period1 = period_fs(i,:);
   period_sec = period(i,:);
   trace2 = delta465_filt(period1(1):period1(2));
   perc90 = prctile(trace2,90);
   perc10 = prctile(trace2,10);
   perc_stitch = perc90-perc10;
   stitch = [stitch perc_stitch];
   trace_2_stitch = [trace_2_stitch trace2];
end

time_stich = 1:1:length(trace_2_stitch );
sec_stich= time_stich/signal_fs;
[pks, pklocs, w, p] = findpeaks(trace_2_stitch, sec_stich, 'MinPeakDistance',10, 'MinPeakProminence',1.5);
figure
findpeaks(trace_2_stitch, sec_stich, 'MinPeakDistance',10, 'MinPeakProminence',1.5);

NE_osc_rate = length(pks)/(length(trace_2_stitch)/signal_fs);
amplitude_mean = mean(stitch);


