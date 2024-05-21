
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
% by Soodeh Moallemian, Ph.D., Rutgers University-2024.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                   STEP-01: Set up the dependencies
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % dep_path is the path to all the dependencies for running this function.
% % You can modify this path based on the path of your system.
% disp('-------------------------------------------------')
% disp('-------------------------------------------------')
% disp('STEP-01: Set up the dependencies')
% disp('-------------------------------------------------')
% disp('-------------------------------------------------')
% disp(['This step is optional, skip this step if you have already added ' ...
%     'the related dependency path.'])
% prompt = ['Have you already added the dependency path? if yes, press "Y" ' ...
%     'and if no, press "N" '];
% Uinput =input(prompt, 's');
% if strcmpi(Uinput,'Y')
%     disp('please enter the path to the dependencies:')
%     new_dep_path = input('Enter dependencies path: ', 's');
%     dep_path = new_dep_path;
% elseif strcmpi(Uinput,'N')
%     disp('please enter the path to the dependencies:')
%     new_dep_path = input('Enter new path: ', 's');
%     dep_path = new_dep_path;
%     addpath(genpath(dep_path));
% end
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   STEP-02: Convert EDF to SET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beep
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
beep
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
beep
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
beep
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-05: Visual quality control of the data')
disp('Here you can check the quality of the signals from each channel.')
disp('you can trim(clip) the data vertically or horizentally.')
disp(['HORIZENTALLY: Please visually check the data and if you find one' ...
    'channel that looks off in more than 50 %% of the cases, ' ...
    'delete that channel by simply clining on the name of the channel ' ...
    'on the left side of the figure. '])
disp(['VERTICALLY: When you visually check the signals and ' ...
'using "1" (the start point of the trim) and "2" the end point'])

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
beep
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
beep
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-07: Bad segment detection (semi-automatic)')
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('********')
disp('This part of the function is semi-automatic.')
disp('********')
fprintf('Please choose the method you want to use for segmenting the bad signal:\n')
method = 'wispic';
N2N3_EEG=csc_artifact_rejection_automated_Eisai(N2N3_EEG, method, 'epoch_length', 6);
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
cd(N2N3_EEG.save_dir);
save(sprintf('%s_N2N3_stagenames_afterfft.mat',sub_id), ...
    'N2N3_stages_afterfft');
%save file with fft
tic
fftfilename = [strrep(newfilename, '_n2n3_noarousals.set', ...
    '_n2n3_noarousals_fftstep'),'.set'];
EEG = pop_saveset(N2N3_EEG,'filename', fftfilename, 'filepath', ...
    sub_fold, 'check', 'on');
disp(['... elapsed time ',num2str(toc/60),' saving set file']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          STEP-08: Check data after artifact rejection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beep
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-08: Check data after artifact rejection')
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('this is your last chance to look into the data and check its quality visually.')
disp(['Note! Now, on the figure, you will only see the N2N3 sleep stages ' ...
    'which represent the nonREM sleep time course. Therefore, it is normal to see ' ...
    'less on the horizental timeline.'])
disp(['if you have not done any triming before, it is best to add two clips' ...
    ' at the end of the signals (for horizental clipping) '])
N2N3_EEG = csc_eeg_plotter_NEW_Eisai(N2N3_EEG);
if ~isempty(N2N3_EEG.hidden_channels)
    N2N3_EEG.bad_channels{1} = cat(2, N2N3_EEG.bad_channels{1}, ...
        str2double({N2N3_EEG.chanlocs(N2N3_EEG.hidden_channels).labels}));
    N2N3_EEG = pop_select(N2N3_EEG, 'nochannel', N2N3_EEG.hidden_channels);
    tic
    N2N3_EEG.bad_channels{1} = cat(2,N2N3_EEG.bad_channels{1},new_bad_channels);
    N2N3_EEG = pop_select(N2N3_EEG, 'nochannel', N2N3_EEG.hidden_channels);
    disp(['... elapsed time ',num2str(toc/60),' removing bad channels']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      STEP-09: Check the Clippings(they must be in even numbers)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beep
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-09: Check the Clippings' )
disp('-------------------------------------------------')
disp('-------------------------------------------------')
clip_check=0;
while clip_check == 0
    if isfield(N2N3_EEG,'csc_event_data')
        if ~isempty(N2N3_EEG.csc_event_data)
            if size(N2N3_EEG.csc_event_data,1) > 2
                for n = 2:size(N2N3_EEG.csc_event_data,1)
                    if isequal(N2N3_EEG.csc_event_data(n,3), N2N3_EEG.csc_event_data(n-1,3))
                        disp('there are back to back same color markings');
                        beep
                        N2N3_EEG = csc_eeg_plotter_NEW_Eisai(N2N3_EEG);
                        clip_check=0;
                    else
                        clip_check=1;
                    end
                end
            end
        end
    else
    end
    check_csc_event_data =  ~mod(size(N2N3_EEG.csc_event_data,1),2);

    if check_csc_event_data == 0
        disp('there is an uneven number of markings');
        beep
        N2N3_EEG = csc_eeg_plotter_NEW_Eisai(N2N3_EEG);
        clip_check=0;
    else
        clip_check=1;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      STEP-10: Remove epochs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beep
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-10: Remove epochs' )
disp('-------------------------------------------------')
disp('-------------------------------------------------')

if clip_check == 1 && ~isempty(N2N3_EEG.csc_event_data)
    event_starts = cellfun(@(x) strcmp(x, 'event 1'), N2N3_EEG.csc_event_data(:, 1));
    
    % sanity check for artifact event markers
    if sum(event_starts) ~= sum(~event_starts)
        fprintf('\nWarning: uneven number of events, check event_data\n');
    end
    
    % use EEGLAB to remove the points
    N2N3_EEG.bad_segments{1} = [cell2mat(N2N3_EEG.csc_event_data(event_starts, 2)), ...
        cell2mat(N2N3_EEG.csc_event_data(~event_starts, 2))];
    
    % convert the timing from seconds to samples
    N2N3_EEG.bad_segments{1} = floor(N2N3_EEG.bad_segments{1} * N2N3_EEG.srate);
    
    % use EEGLAB to remove the regions
    N2N3_EEG = pop_select(N2N3_EEG, 'nopoint', N2N3_EEG.bad_segments{1});
    
    
end
% save the file with channels/epochs deleted & marked in n2n3 file
tic
n2n3filename = [strrep(fftfilename,'_n2n3_noarousals_fftstep.set','_n2n3_noarousals_fftstep_artifcorr'),'.set'];
N2N3_EEG = pop_saveset(N2N3_EEG,'filename', n2n3filename, 'filepath', N2N3_EEG.save_dir, 'check', 'on');

load(sprintf('%s_N2N3_stagenames_afterfft',sub_id));
numsamps_afterfft = zeros(1,length(N2N3_stages_afterfft));
bad_segments = cell2mat(N2N3_EEG.bad_segments);

for j = 1:size(bad_segments,1)
    numsamps_afterfft(bad_segments(j,1):bad_segments(j,2)) = 1;
end

N2N3_stages_afterfft_snipped = N2N3_stages_afterfft(numsamps_afterfft ==0);
% Removing the last sample, since this is a mismatch, can recheck this
% later
if size(N2N3_stages_afterfft_snipped,2) ~= size(N2N3_EEG.data,2)
    N2N3_stages_afterfft_snipped(end) = [];
end
N2N3_EEG.N2N3_stagenames_afterfft_snipped =  N2N3_stages_afterfft_snipped;
save([N2N3_EEG.save_dir sub_id '_N2N3_stagenames_afterfft_snipped.mat'],'N2N3_stages_afterfft_snipped');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          STEP-12: Check bad channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-12: Check bad channels')
disp('-------------------------------------------------')
disp('-------------------------------------------------')
finished = 'N';
while ~isempty(N2N3_EEG.hidden_channels) || strcmp(finished,'N')
    N2N3_EEG = csc_eeg_plotter_NEW_Eisai(N2N3_EEG);
       
    new_bad_channels = str2double({N2N3_EEG.chanlocs(N2N3_EEG.hidden_channels).labels});
    
    tic
    N2N3_EEG.bad_channels{1} = cat(2,N2N3_EEG.bad_channels{1},new_bad_channels);
    if ~isempty(N2N3_EEG.hidden_channels)
        N2N3_EEG = pop_select(N2N3_EEG, 'nochannel', N2N3_EEG.hidden_channels);
        disp(['... elapsed time ',num2str(toc/60),' removing bad channels']);
%       N2N3_EEG = pop_reref( N2N3_EEG, [] );
        N2N3_EEG.hidden_channels = [];
    end
   if isfield(N2N3_EEG,'psd'); N2N3_EEG = rmfield(N2N3_EEG, 'psd');end
    for n = 1:N2N3_EEG.trials
        [N2N3_EEG.psd.data(:,:,:,n),N2N3_EEG.psd.Hzbins] = psddata(squeeze(N2N3_EEG.data(:,:,n)),N2N3_EEG.srate,6,40,1); % only use fast for massive data
      
    end
    close all;
    plot_bands_spectra_EEG_PAMMS_Eisai(N2N3_EEG);
    corr_ans = 0;
    while ~corr_ans
        finished = input('Finished Y/N? [Y]','s');
        if ~ismember(finished, {'Y','N'})
            corr_ans = 0;
        else
            corr_ans = 1;
        end
    end
end

%save file data set
tic
if finished == 'Y'
    finishedfilename = [strrep(n2n3filename,'_n2n3_noarousals_fftstep_artifcorr.set','_n2n3_noarousals_fftstep_artifcorr_final'),'.set'];
    %           finishedfilename = [strrep(EEG.filename,'_n2n3_noarousals_artifcorr_fftstep.set','_n2n3_noarousals_artifcorr_fftstep_final'),'.set'];
    N2N3_EEG.datfile = finishedfilename;
    N2N3_EEG = pop_saveset(N2N3_EEG,'filename', finishedfilename, 'filepath', sub_fold, 'check', 'on');
end
disp(['... elapsed time ',num2str(toc/60),' saving set file']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          STEP-13: Segmenting the N2 and N3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beep
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-13: Segmenting the N2 and N3')
disp('-------------------------------------------------')
disp('-------------------------------------------------')
% for the moment we are not rejecting any channels. so in this function we
% don't perform the removing bad channels part. We will only save the final
% data after checking the atrifact rejection from N2 and N3 seperately

load(fullfile(sub_fold,[sub_id '_N2N3_stagenames_afterfft.mat'])); 
N2_EEG = N2N3_EEG;

finishedfilename = 'preprocessed_eeg';
N2cleaned_filename = strrep(finishedfilename,'.set','_N2only.set');
N2only_cleaned = (strcmp(N2N3_stages_afterfft,'N2'));
swa_selectStagesEEGLAB_BAR(N2_EEG, N2only_cleaned,'N2only_cleaned',N2cleaned_filename);
N2only_cleaned = pop_saveset(N2N3_EEG, 'filename', N2cleaned_filename, 'filepath', sub_fold ,'check', 'on');
 

N3_EEG = N2N3_EEG;
N3cleaned_filename = strrep(finishedfilename,'.set','_N3only.set');

N3only_cleaned = (strcmp(N2N3_stages_afterfft,'N3'));
swa_selectStagesEEGLAB_BAR(N3_EEG, N3only_cleaned,'N3only_cleaned',N3cleaned_filename);
N3only_cleaned = pop_saveset(N2N3_EEG, 'filename', N3cleaned_filename, 'filepath',sub_fold ,'check', 'on');



NREMpsd = N2N3_EEG.psd;
save NREMpsd NREMpsd;

end