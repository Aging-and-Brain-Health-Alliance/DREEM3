% %----------------------------------------------------------------------%
% name:          run_getModulationIndexLocal
%
%       Script for calculating phase-amplitude coupling measures for Eisai 
%       between two frequency bands of interest: Mean Vector Length (what Canolty
%       refers to as Modulation Index, but this conflicts with other
%       literature and MVL is a natural name for this, since it is the mean
%       absolute length of a complex vector); and Modulation Index
%       (introduced by Tort, and based on the statistics of the histogram
%       of amplitudes within phase bins, specifically a normalized
%       Kullback-Leibler divergence between observed distribution and a
%       uniform one) - here MI is modified to cover 2 cycles! in order to
%       distinguish events pre-trough vs post-trough.
%
%       This version computes the phase-amplitude coupling measures only
%       for the signal cut-outs around Slow Oscillations (SOs) and not for
%       a continuous signal across the whole stage. This was motivated by
%       the fact that the interest was to examine proper slow oscillations,
%       while taking into account the rest of the signal might create
%       spurious coupling, or obscure the existing one (in cases when SOs
%       are not a majority of the signal). 
%
%       First, SOs are detected with a simple threshold-based algorithm.
%       EEG is filtered into 0.1-4 Hz, consecutive zero crossings are
%       identified, and then deflections of sufficient depth and width are
%       marked as slow oscillations. The filtering is in-built, since most
%       other studies detect SOs on filtered EEG traces. Detection criteria
%       are customizable by the user; suggested parameters are here similar
%       to Massiminis paper. 
%
%       Then filtered and hilbert-transformed phase-providing and
%       amplitude-providing signals are passed to the Modulation Index
%       function along with the sample numbers of detected SOs. Cut-outs
%       around the trough are chosen based on the cycle length, so it is
%       locally adaptive to the SO shape, rather than a pre-defined window
%       (which might result in over-representation of some phase bins,
%       technically): one full cycle before the trough, one full cycle
%       after the trough. The function computes raw values of MVL and MI,
%       and also performs surrogate data generation, so that normalized
%       values of MVL and MI can be obtained, as well as p-values. For
%       maximum customisability, it writes out the surrogate distributions
%       of MVL and MI, as well as computed / estimated angles of the
%       coupling and their distributions.
%       For more detail, see description of getModulationIndexLocal.
%
%       Script saves all the parameters, time of computation, current date
%       and time at the finish, all phase-amp coupling measures and
%       relevant distributions, as well as detected slow oscillations
%
% dependency:   - format of staging file
%               - format of raw data file
%               - naming of the files 
%               - function: edfread
%               - function: detectSO
%               - function: chooseEventsBySleepStage
%               - function: getModulationIndexLocal
%
% WARNING: current version is aimed at nap non-dense EEG data, because it
%          reads in the whole record at once, and with too long recordings
%          or too many electrodes it will likely run out of memory; for
%          other data has to be updated to work via channel-by-channel
%          read-in
%
%          current version is suited for SLOW OSCILLATIONS as
%          phase-providing band ONLY, because the function
%          getModulationIndexLocal assumes the same frequency for for
%          cutting out 2 cycles around SO trough, and for phase-providing;
%          for other phase-providing bands has to be updated to have an
%          additional input frequency for example
%
%           
% %----------------------------------------------------------------------%


% Here are modifiable parameters:
% ----------------------------------------------------------------------- %
clc
clear
close all

start_time = datetime("now")


addpath(genpath('C:\Users\mt2294\OneDrive - University of Bath\Documents\MATLAB\eeglab2023.0'));
addpath(genpath('C:\Users\mt2294\OneDrive - University of Bath\UCI\Data analysis scripts\UWMpipeline\Scripts\Post processing'));
eeglab

plotflag=0;

Data_dir='C:\Users\mt2294\OneDrive - University of Bath\UCI\Data analysis scripts\UWMpipeline\Scripts\RESTED\temp_loop\';
Subjects={'RESTED108_D1'};

% filter frequencies used:
filtfreq_slow = [0.5 1.25];                % lower and upper frequency bound [Hz] for the phase-providing band
filtfreq_fast = [13 16];                  % lower and upper frequency bound [Hz] for the amplitude-providing band

% size of the epochs used for staging [sec]:
epoch_size = 30;

% SO detection criteria:
so_volt_threshold = 80;          % minimum allowable depth of the negative half of the SO [uV]
so_peak2peak = 80;               % minimum allowable distance from deepest to highest voltage of the SO [uV]
so_width_min = 0.3;              % minimum allowable width of the negative half of the SO [sec] (distance first and second zero crossing)
so_width_max = 1.5;                % maximum allowable width of the negative half of the SO [sec]
so_width_max_2 = 10;             % maximum allowable width of the SO with peak [sec] (distance between first and third zero crossing)

% filter orders to use for butterworth filter design, these have been
% manually adjusted for the SO and sigma bands:
fo_slow_high = 4;
fo_slow_low = 4;
fo_fast_high = 9;   %9;
fo_fast_low = 8;   %8;

% number of times a surrogate dataset is constructed:
N_shift = 2000;

for sub=1:length(Subjects)

    EEG = pop_loadset('filename',[ Subjects{sub}  '_SnippedSleep_filtered_bandpass_n2n3_noarousals_fftstep_artifcorr_final.set'],'filepath',[Data_dir Subjects{sub} '\']);

    n_channel=EEG.nbchan;


    % prepare space for output variables:
    events_slow = cell(n_channel,1);
    mod_c_raw = NaN*ones(n_channel,1);  % mean vector length, raw
    mod_t_raw = NaN*ones(n_channel,1);  % modulation index, raw
    mod_c_norm = NaN*ones(n_channel,1);    % mean vector length, normalized
    mod_t_norm = NaN*ones(n_channel,1);    % modulation index, normalized
    mod_c_pval = NaN*ones(n_channel,1);   % probability of obtaining such a high MVL as compared to surrogate
    mod_t_pval = NaN*ones(n_channel,1);   % probability of obtaining such a high MI as compared to surrogate
    mod_c_distr = cell(n_channel,1);    % distribution of surrogate values of MVL
    mod_t_distr = cell(n_channel,1);    % distribution of surrogate values of MI
    mod_c_angle = NaN*ones(n_channel,1);   % the angle of the complex vector, i.e. phase at which amplitude tends to peak
    mod_t_hist = cell(n_channel,1);      % histogram of amplitudes of fast signal within bins -2pi:pi/18:2pi, centered on SO trough


    electrodes={};

    for ch=1:n_channel
        electrodes{ch}=EEG.chanlocs(ch).labels;
        [slow_events,~] = detectSO(double(EEG.data(ch,:)),EEG.srate,so_volt_threshold,so_peak2peak,so_width_min,so_width_max,so_width_max_2);
        if ~isempty(slow_events)
            slow_events = slow_events(:,2);
        end
        events_slow{ch} = slow_events;
        % slow band:
        [b1_slow,a1_slow] = butter(fo_slow_high,filtfreq_slow(1)./(EEG.srate/2),'high');
        [b2_slow,a2_slow] = butter(fo_slow_low,filtfreq_slow(2)./(EEG.srate/2),'low');
        x = filtfilt(b1_slow,a1_slow,double(EEG.data(ch,:)));
        channel_slow = filtfilt(b2_slow,a2_slow,x);
        % fast band:
        [b1_fast,a1_fast] = butter(fo_fast_high,filtfreq_fast(1)./(EEG.srate/2),'high');
        [b2_fast,a2_fast] = butter(fo_fast_low,filtfreq_fast(2)./(EEG.srate/2),'low');
        x = filtfilt(b1_fast,a1_fast,double(EEG.data(ch,:)));
        channel_fast = filtfilt(b2_fast,a2_fast,x);

        [mi_raw_c,mi_matrix_c,mi_angle,mi_raw_t,mi_matrix_t,mi_hist,edges] = getModulationIndexLocal_NS(hilbert(channel_slow),hilbert(channel_fast),slow_events,EEG.srate,N_shift,plotflag);

        mod_c_raw(ch) = mi_raw_c;
        mod_c_distr{ch} = mi_matrix_c;
        mod_c_angle(ch) = mi_angle;
        mod_c_pval(ch) = length(find(mi_matrix_c > mi_raw_c))./length(mi_matrix_c);
        [mm,ss] = normfit(mi_matrix_c);
        mod_c_norm(ch) = (mi_raw_c - mm)./ss;

        mod_t_raw(ch) = mi_raw_t;
        mod_t_distr{ch} = mi_matrix_t;
        mod_t_hist{ch} = mi_hist;
        mod_t_pval(ch) = length(find(mi_matrix_t > mi_raw_t))./length(mi_matrix_t);
        [mm,ss] = normfit(mi_matrix_t);
        mod_t_norm(ch) = (mi_raw_t - mm)./ss;




    end

    file_out=[Data_dir Subjects{sub} '\' Subjects{sub} '_ModulationIndexRe'];

    save(file_out,...
        'electrodes',...
        'so_volt_threshold',...
        'so_peak2peak',...
        'so_width_min',...
        'so_width_max',...
        'so_width_max_2',...
        'events_slow',...
        'filtfreq_slow',...
        'filtfreq_fast',...
        'events_slow',...
        'N_shift',...
        'mod_c_raw',...
        'mod_t_raw',...
        'mod_c_norm',...
        'mod_t_norm',...
        'mod_c_distr',...
        'mod_t_distr',...
        'mod_c_angle',...
        'mod_t_hist',...
        'mod_c_pval',...
        'mod_t_pval',...
        'edges')
end

