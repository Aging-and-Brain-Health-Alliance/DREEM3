function main_eeg_preprocess(sub_id, source_path)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function prepares EEG raw data, performs preprocessing, and saves
% files needed for further investigation (e.g., sleep spindles detection).
% To run this function, eeglab must be installed. The function will save
% the generated files in the subject's folder 
%
% INPUTS:
%   sub_id: Subject ID
%   source_path: Path to the source folder containing EEG data
%
% EXAMPLE INPUTS:
%   sub_id = 'COV207';
%   source_path = '/Users/smoallemian/Desktop/Dreem_study/data/';
%   dest_path = source_path
%
% You can modify this function's paths for your installation.
%
% This code was written by Soodeh Moallemian, Ph.D., Rutgers University-2024.
% For questions, contact s.moallemian@rutgers.edu.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   STEP-01: Set up the dependencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dep_path is the path to all the dependencies for running this function.
% You can modify this path based on the path of your system.
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-01: Set up the dependencies')
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('*********')
disp('Is this the path you have saved the dependencies?')
dep_path = '/Users/smoallemian/Desktop/Dreem_study/Dreem_Scripts_NS';
fprintf('%s\n', dep_path)
disp('*********')

prompt = 'If YES, please press enter(return) button. Otherwise press "N": ';
disp('*********')
user_input = input(prompt, 's');

if strcmpi(user_input, 'N')
    fprintf('Please enter your desired path:\nExample: /Users/XXXXX/Desktop/XXXX\n')
    new_dep_path = input('Enter new path: ', 's');
    addpath(genpath(new_dep_path))
else
    addpath(genpath(dep_path))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   STEP-02: Convert EDF to SET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-02: Convert EDF to SET')
disp('-------------------------------------------------')
disp('-------------------------------------------------')
sub_fold = fullfile(source_path, sub_id);
cd(sub_fold)
eeglab

%read the edf file and load it to eeglab
edf_file = spm_select('FPList', sub_fold, '.*\.edf');
edf_file = pop_biosig(edf_file);

% get the frequency rate for future calculations
Fs = edf_file.srate;

%read epochs from the txt report file in a table
report_file = spm_select('FPList', sub_fold, '.*\.txt');
report_tab = readtable(report_file);

% Read the duration of each epoch which is 30 seconds
% Define the epochs based on the sleep stages in txt report

epoch_duration = report_tab.Duration_s_; 

% Mapping sleep stages
% Define the sleep stages and their related score
%  0 == wake(sleep-s0), 1/2/3 == NREM, 4 == REM, 5 == MT
sleep_stages = containers.Map({'SLEEP-S0', 'SLEEP-S1', 'SLEEP-S2', 'SLEEP-S3', 'SLEEP-REM', 'SLEEP-MT'}, {0, 1, 2, 3, 4, 5});

% Read the epoch based on the report file and sleep stages
epoch = cellfun(@(x) sleep_stages(x), report_tab.SleepStage);
epoch = epoch';
% We add 1 to each epoch to be able to find the first 0 because "find" 
% function searches for non-zero. SO FROM NOW ON< WAKE IS 1.
epoch = epoch+1;
% Find the first and last epoches
first_epoch = find(epoch == 1, 1, 'first');
last_epoch = find(epoch == 1, 1, 'last');
sleep_stages = epoch(first_epoch:last_epoch);

% Define a name for the raw file and create a .set file (readable for eeglab)
filename = [sub_id 'RawSleep.set'];
% Define Channel names for DREEM3
EEG = pop_select(edf_file, 'channel', {'F7-O1', 'F8-O2', 'F8-F7', 'F8-O1', 'F7-O2'});
% Save the updated file
EEG = pop_saveset(EEG, 'filename', filename, 'filepath', sub_fold, 'check', 'on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             STEP-03: Clip EDF based on sleep stages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-03: Clip EDF based on sleep stages')
disp('-------------------------------------------------')
disp('-------------------------------------------------')
% First data is calculated based on the frequency rate (FS)
% * the duration of each epoch(which is 30) using the real first epoch 
% (before adding the 1)
first_data = (first_epoch - 1) * Fs * epoch_duration(1) + 1;
%first data is calculated based on the frequency rate (FS)
% * 30(which is the duration of each epoch) using the last epoch 
last_data = last_epoch * Fs * epoch_duration(1);
EEG = pop_select(EEG, 'point', [first_data last_data]);
EEG.sleepstage_labels = sleep_stages;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          STEP-04: Map sleep stages to all samples of sleep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-04: Map sleep stages to all samples of sleep')
disp('-------------------------------------------------')
disp('-------------------------------------------------')

sleep_stages2 = repmat(sleep_stages, 1, Fs * 30);
EEG.stagenames = reshape(sleep_stages2',size(sleep_stages2,1)*size(sleep_stages2,2),1)-1;
stagenames = containers.Map({0, 1, 2, 3, 4, 5},{'W', 'N1', 'N2', 'N3','REM', 'MT'});
temp = cellfun(@(x) stagenames(x), num2cell(EEG.stagenames), 'UniformOutput', false);
EEG.stagenames = temp;

filename = sprintf('%s_SnippedSleep_filtered_bandpass.set', sub_id);
EEG = pop_saveset(EEG, 'filename', filename, 'filepath', sub_fold, 'check', 'on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          STEP-05: Visual quality control of the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-05: Visual quality control of the data')
disp('-------------------------------------------------')
disp('-------------------------------------------------')
EEG = pop_loadset(fullfile(sub_fold,[sub_id '_SnippedSleep_filtered_bandpass.set']));
load(fullfile (dep_path , 'Dependencies', 'chanlocs_dreem3_mock.mat'), 'chanlocs1');
EEG.chanlocs = chanlocs1;
EEG.original_chanlocs = EEG.chanlocs;
EEG = csc_eeg_plotter_NEW_Eisai(EEG);
EEG.bad_channels{1} = EEG.hidden_channels;
EEG = pop_select(EEG, 'nochannel', EEG.bad_channels{1});
% save original bc data
EEG.save_dir= sub_fold;
EEG.filename = strrep(EEG.filename,'.set','_bc1pass.set');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 STEP-06: Create N2N3 set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-06: Create N2N3 set')
disp('-------------------------------------------------')
disp('-------------------------------------------------')

N2N3_samples_inds = [];
N2only_inds=[];
N3only_inds=[];

for i =1:length(EEG.stagenames)
    if strcmp(EEG.stagenames{i},'N2')||strcmp(EEG.stagenames{i},'N3')
        N2N3_samples_inds(i)=1;
    end
    if strcmp(EEG.stagenames{i},'N2')
        N2only_inds(i)=1;
    end
    if strcmp(EEG.stagenames{i},'N3')
        N3only_inds(i)=1;
    end
end
t_n2n3=find(N2N3_samples_inds); % find the indices both for N2 and N3
dif=find(t_n2n3(2:end)-t_n2n3(1:end-1)>1);
ends=[t_n2n3(dif),t_n2n3(end)];
starts=[t_n2n3(1),t_n2n3(dif+1)];
n2n3_intervals=[starts',ends'];

N2N3_EEG=pop_select( EEG, 'point',n2n3_intervals );

newfilename  = [strrep(N2N3_EEG.filename,'_bc1pass.set','_n2n3_noarousals'),'.set'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       STEP-07: Bad segment detection (semi-automatic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-07: Bad segment detection (semi-automatic)')
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('********')
disp('This part of the function is semi-automatic.')
disp('********')
fprintf('Please choose the method you want to use for segmenting the bad signal:\n')
prompt = 'Press "w" for wispic or "e" for eeglab\n';
user_input = input(prompt , 's');
if strcmpi(user_input, 'w')
   method = 'wispic';
    N2N3_EEG=csc_artifact_rejection_automated_Eisai(N2N3_EEG, method, 'epoch_length', 30);
    N2N3_EEG = pop_select(N2N3_EEG, 'notime', N2N3_EEG.bad_regions);
    N2N3_stages=N2N3_EEG.stagenames(find(N2N3_samples_inds));
    N2N3_artifactsamps = zeros(1,length(N2N3_stages));
    bad_regions = N2N3_EEG.bad_regions;

    for i = 1:length(bad_regions)
        start = bad_regions(i,1)*N2N3_EEG.srate;
        last = bad_regions(i,2)*N2N3_EEG.srate;
        if start == 0
            start = 1;
        end
        N2N3_artifactsamps(start:last) = 1;
    end

    N2N3_stages_afterfft = N2N3_stages(N2N3_artifactsamps ==0);
    save(sprintf('%s_N2N3_stagenames_afterfft.mat',sub_id),'N2N3_stages_afterfft');
%save file with fft
    tic
    fftfilename = [sprintf([newfilename '_n2n3_noarousals.set','_n2n3_noarousals_fftstep'],'.set')];
    EEG = pop_saveset(N2N3_EEG,'filename', fftfilename, 'filepath', sub_fold, 'check', 'on');
    disp(['... elapsed time ',num2str(toc/60),' saving set file']);
else

    fprintf(['THIS PART NEEDS MORE WORK NOW. Please contact us for the updates. ' ...
        'at the moment the function works only using wispic method'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          STEP-08: Check data after artifact rejection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-08: Check data after artifact rejection')
disp('-------------------------------------------------')
disp('-------------------------------------------------')
N2N3_EEG = csc_eeg_plotter_NEW_Eisai(N2N3_EEG);
if ~isempty(N2N3_EEG.hidden_channels)
    N2N3_EEG.bad_channels{1} = cat(2, N2N3_EEG.bad_channels{1}, str2double({N2N3_EEG.chanlocs(N2N3_EEG.hidden_channels).labels}));
    N2N3_EEG = pop_select(N2N3_EEG, 'nochannel', N2N3_EEG.hidden_channels);
    tic
    N2N3_EEG.bad_channels{1} = cat(2,N2N3_EEG.bad_channels{1},new_bad_channels);
    N2N3_EEG = pop_select(N2N3_EEG, 'nochannel', N2N3_EEG.hidden_channels);
    disp(['... elapsed time ',num2str(toc/60),' removing bad channels']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          STEP-09: Segmenting the N2 and N3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-09: Segmenting the N2 and N3')
disp('-------------------------------------------------')
disp('-------------------------------------------------')
% for the moment we are not rejecting any channels. so in this function we
% don't perform the removing bad channels part. We will only save the final
% data after checking the atrifact rejection from N2 and N3 seperately

load(fullfile(sub_fold,[sub_id '_N2N3_stagenames_afterfft.mat'])); 
N2_EEG = N2N3_EEG;

finishedfilename = 'preprocessed_eeg'
N2cleaned_filename = strrep(finishedfilename,'.set','_N2only.set');
N2only_cleaned = (strcmp(N2N3_stages_afterfft,'N2'));
swa_selectStagesEEGLAB_BAR(N2_EEG, N2only_cleaned,'N2only_cleaned',N2cleaned_filename);
N2only_cleaned = pop_saveset(N2N3_EEG, 'filename', N2cleaned_filename, 'check', 'on');
 

N3_EEG = N2N3_EEG;
N3cleaned_filename = strrep(finishedfilename,'.set','_N3only.set');

N3only_cleaned = (strcmp(N2N3_stages_afterfft,'N3'));
swa_selectStagesEEGLAB_BAR(N3_EEG, N3only_cleaned,'N3only_cleaned',N3cleaned_filename);
N3only_cleaned = pop_saveset(N2N3_EEG, 'filename', N3cleaned_filename, 'check', 'on');



end
