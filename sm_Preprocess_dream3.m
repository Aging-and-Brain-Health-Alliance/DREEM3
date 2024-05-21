close all
clear 
clc
addpath('/Users/smoallemian/Desktop/Dreem_study/Dreem_Scripts_NS');
eeglab;
BasePath = fullfile('/Users/smoallemian/Desktop/Dreem_study/data/');
cd (BasePath)
% Define the subject count (assuming you have 10 subjects, for example)
sub_id = 'COV207'; % Change this value based on the subject you want to process
finished = 'N';
%% edf2set conversion
sub_fold = fullfile(BasePath,sub_id);
cd (sub_fold)

%read the edf file and load it to eeg lab
edf_file = spm_select('FPList',sub_fold, '.*.edf');
EEG = pop_biosig(edf_file);
%the frequency rate
Fs = EEG.srate; %250

%read epochs from the txt report file in a table

report_file = spm_select('FPList',sub_fold, '.*.txt');
report_tab = readtable(report_file);
epoch_duration = report_tab.Duration_s_; % read the duration of each epoch which is 30 seconds
%define the epochs based on the sleep stages in txt report
%  0 == wake, 1/2/3 == NREM, 4 == REM, 5 == MT

sleep_stages = {'SLEEP-S0', 'SLEEP-S1', 'SLEEP-S2', 'LEEP-S3',...
    'SLEEP-REM', 'SLEEP-MT'};

epoch = [];
for ii = 1: size(report_tab,1)
    SleepStage = report_tab.SleepStage{ii};
switch SleepStage 
    case sleep_stages{1}
        epoch(end+1) = 0; %  0 == wake
    case sleep_stages{2}
        epoch(end+1) = 1; %  1 == nonREM
    case sleep_stages{3}
        epoch(end+1) = 2;%  2 == nonREM
    case sleep_stages{4}
        epoch(end+1) = 3;%  3 == nonREM
    case sleep_stages{5}
        epoch(end+1) = 4;%  4 == REM
    case sleep_stages{6}
    epoch(end+1) = 5;%  5 == MT
end
end
% add 1 to all epochs to skip the wake stage
epoch = epoch+1;

%finding the first epoch
first_epoch= find(epoch==1,1,'first');%first wake epch that now is 1 because we added 1
%finding the last ep[och
last_epoch= find (epoch==1,1,'last');%last wake epoch

sleep_stages=epoch(first_epoch:last_epoch);%in between first to last stage of scored sleep


%% save and creat .set
filename = [sub_id '_RawSleep.set'];

%define the channel names
CHs={'F7-O1','F8-O2','F8-F7','F8-O1','F7-O2'};
CH_inds=[];
CHs_name={};


for ch=1:length(CHs)
    for ii = 1:length(EEG.chanlocs)
        if  strcmp(EEG.chanlocs(ii).labels, CHs{ch})
            CH_inds(ch)=ii;

            CHs_name{ch}=EEG.chanlocs(ii).labels;
        end
    end
end

EEG = pop_select( EEG, 'channel',CHs_name);

EEG = pop_saveset(EEG, 'filename', filename, 'filepath', sub_fold, 'check', 'on');

%% clip the EDF based on sleep stages

first_data=(first_epoch-1)*Fs*30+1;
last_data=last_epoch*Fs*30;
EEG = pop_select( EEG, 'point',[first_data last_data] );
EEG.sleepstage_lables=sleep_stages;

filenameNew='snippedData';
%% mapping sleep stages to all samples of sleep


sleep_stages2=repmat(sleep_stages,1,Fs*30);
sleep_stages2=reshape(sleep_stages2',size(sleep_stages2,1)*size(sleep_stages2,2),1)-1;
EEG.stagenames = sleep_stages2;
theCells = cell(size(EEG.stagenames));
theCells(EEG.stagenames==0) = {'W'};
theCells(EEG.stagenames==1) = {'N1'};
theCells(EEG.stagenames==2) = {'N2'};
theCells(EEG.stagenames==3) = {'N3'};
theCells(EEG.stagenames==5) = {'REM'};
theCells(EEG.stagenames==6) = {'MT'};
%theCells(EEG.stagenames==7) = {'NS'};
EEG.stagenames = theCells;

filename = sprintf('%s_SnippedSleep_filtered_bandpass.set',sub_id);
EEG = pop_saveset(EEG, 'filename', filename, 'filepath', sub_fold, 'check', 'on');


%% Data quality control (to check if there are some very bad channels missing, check visually)
EEG = pop_loadset([sub_id '_SnippedSleep_filtered_bandpass.set']);
addpath '/Users/smoallemian/Desktop/Dreem_study/Dreem_Scripts_NS/Dependencies/csc-eeg-tools-develop'
load('/Users/smoallemian/Desktop/Dreem_study/Dreem_Scripts_NS/Dependencies/chanlocs_dreem3_mock.mat','chanlocs1')
EEG.chanlocs=chanlocs1;
EEG.original_chanlocs = EEG.chanlocs;
EEG = csc_eeg_plotter_NEW_Eisai(EEG);
EEG.bad_channels{1} = EEG.hidden_channels;
EEG = pop_select(EEG, 'nochannel', EEG.bad_channels{1});
% save original bc data
EEG.save_dir= sub_fold;
EEG.filename = strrep(EEG.filename,'.set','_bc1pass.set');
%EEG = pop_saveset(EEG,'filename', EEG.filename, 'filepath', EEG.save_dir); %'check', 'on'


%% create N2N3 set

%load(sprintf('%s_EEGevents_corrected.mat',subj_id));
N2N3_samples_inds=[];
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
t_n2n3=find(N2N3_samples_inds); % find the indexes both for N2 and N3
dif=find(t_n2n3(2:end)-t_n2n3(1:end-1)>1);
ends=[t_n2n3(dif),t_n2n3(end)];
starts=[t_n2n3(1),t_n2n3(dif+1)];
n2n3_intervals=[starts',ends'];

N2N3_EEG=pop_select( EEG, 'point',n2n3_intervals );

newfilename  = [strrep(N2N3_EEG.filename,'_bc1pass.set','_n2n3_noarousals'),'.set'];

%% semi-automatic bad segment detection
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% change the threshold to higher or for automatic just accept
addpath '/Users/smoallemian/Desktop/Dreem_study/Dreem_Scripts_NS/Dependencies/csc-eeg-tools-develop/preprocessing/artifact_rejection'
addpath '/Users/smoallemian/Desktop/Dreem_study/Dreem_Scripts_NS/Dependencies/csc-eeg-tools-develop/dependencies'
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
cd(N2N3_EEG.save_dir); %made this change on 12/29 so the file would be save in temp instead of CRSP (by AD)
save(sprintf('%s_N2N3_stagenames_afterfft.mat',sub_id),'N2N3_stages_afterfft');
%save file with fft
tic
fftfilename = [strrep(newfilename,'_n2n3_noarousals.set','_n2n3_noarousals_fftstep'),'.set'];
EEG = pop_saveset(N2N3_EEG,'filename', fftfilename, 'filepath', EEG.save_dir, 'check', 'on');
disp(['... elapsed time ',num2str(toc/60),' saving set file']);



  %% check data after artifact rejection - if there are still bad segments, id
% or bad channels id them
%make sure to thoroughly check for anymore bad segments, as you will
%not be able to snip out anymore
N2N3_EEG = csc_eeg_plotter_NEW_Eisai(N2N3_EEG);
if ~isempty(N2N3_EEG.hidden_channels)
    new_bad_channels = str2double({N2N3_EEG.chanlocs(N2N3_EEG.hidden_channels).labels});
    
    tic
    N2N3_EEG.bad_channels{1} = cat(2,N2N3_EEG.bad_channels{1},new_bad_channels);
    N2N3_EEG = pop_select(N2N3_EEG, 'nochannel', N2N3_EEG.hidden_channels);
    disp(['... elapsed time ',num2str(toc/60),' removing bad channels']);
end

%% check that there is an even amount of clippings
%%%%%I have commented this part as we don't have any csc_event_data
% clip_check=0;
% while clip_check == 0
%     if isfield(N2N3_EEG,'csc_event_data')
%         if ~isempty(N2N3_EEG.csc_event_data)
%             if size(N2N3_EEG.csc_event_data,1) > 2
%                 for n = 2:size(N2N3_EEG.csc_event_data,1)
%                     if isequal(N2N3_EEG.csc_event_data(n,3), N2N3_EEG.csc_event_data(n-1,3))
%                         disp('there are back to back same color markings');
%                         beep
%                         N2N3_EEG = csc_eeg_plotter_NEW_Eisai(N2N3_EEG);
%                         clip_check=0;
%                     else
%                         clip_check=1;
%                     end
%                 end
%             end
%         end
%     end
%     check_csc_event_data =  ~mod(size(N2N3_EEG.csc_event_data,1),2);
%     
%     if check_csc_event_data == 0
%         disp('there is an uneven number of markings');
%         beep
%         N2N3_EEG = csc_eeg_plotter_NEW_Eisai(N2N3_EEG);
%         clip_check=0;
%     else
%         clip_check=1;
%     end
%     
% end
% sprintf('run next section')
%% remove epochs
%we don't have badsegments
% % % if clip_check == 1 && ~isempty(N2N3_EEG.csc_event_data)
% % %     event_starts = cellfun(@(x) strcmp(x, 'event 1'), N2N3_EEG.csc_event_data(:, 1));
% % %     
% % %     % sanity check for artifact event markers
% % %     if sum(event_starts) ~= sum(~event_starts)
% % %         fprintf('\nWarning: uneven number of events, check event_data\n');
% % %     end
% % %     
% % %     % use EEGLAB to remove the points
% % %     N2N3_EEG.bad_segments{1} = [cell2mat(N2N3_EEG.csc_event_data(event_starts, 2)), ...
% % %         cell2mat(N2N3_EEG.csc_event_data(~event_starts, 2))];
% % %     
% % %     % convert the timing from seconds to samples
% % %     N2N3_EEG.bad_segments{1} = floor(N2N3_EEG.bad_segments{1} * N2N3_EEG.srate);
% % %     
% % %     % use EEGLAB to remove the regions
% % %     N2N3_EEG = pop_select(N2N3_EEG, 'nopoint', N2N3_EEG.bad_segments{1});
% % %     %
% % %     
% % % end
% % % % save the file with channels/epochs deleted & marked in n2n3 file
% % % tic
% % % n2n3filename = [strrep(fftfilename,'_n2n3_noarousals_fftstep.set','_n2n3_noarousals_fftstep_artifcorr'),'.set'];
% % % N2N3_EEG = pop_saveset(N2N3_EEG,'filename', n2n3filename, 'filepath', N2N3_EEG.save_dir, 'check', 'on');
% % % load(sprintf('%s_N2N3_stagenames_afterfft',sub_id));
% % % numsamps_afterfft = zeros(1,length(N2N3_stages_afterfft));
% % % bad_segments = cell2mat(N2N3_EEG.bad_segments);
% % % 
% % % for j = 1:size(bad_segments,1)
% % %     numsamps_afterfft(bad_segments(j,1):bad_segments(j,2)) = 1;
% % % end

N2N3_stages_afterfft_snipped = N2N3_stages_afterfft(numsamps_afterfft ==0);
% Removing the last sample, since this is a mismatch, can recheck this
% later
if size(N2N3_stages_afterfft_snipped,2) ~= size(N2N3_EEG.data,2)
    N2N3_stages_afterfft_snipped(end) = [];
end
N2N3_EEG.N2N3_stagenames_afterfft_snipped =  N2N3_stages_afterfft_snipped;
save([N2N3_EEG.save_dir sub_id '_N2N3_stagenames_afterfft_snipped.mat'],'N2N3_stages_afterfft_snipped');
sprintf('run next section')
          %% this while loop will not take any bad segments out only bad channels;
while ~isempty(N2N3_EEG.hidden_channels) || strcmp(finished,'N')
    N2N3_EEG = csc_eeg_plotter_NEW_Eisai(N2N3_EEG);
       
%     new_bad_channels = str2double({N2N3_EEG.chanlocs(N2N3_EEG.hidden_channels).labels});
%     
%     tic
%     N2N3_EEG.bad_channels{1} = cat(2,N2N3_EEG.bad_channels{1},new_bad_channels);
%     if ~isempty(N2N3_EEG.hidden_channels)
%         N2N3_EEG = pop_select(N2N3_EEG, 'nochannel', N2N3_EEG.hidden_channels);
%         disp(['... elapsed time ',num2str(toc/60),' removing bad channels']);
% %         N2N3_EEG = pop_reref( N2N3_EEG, [] );
%         N2N3_EEG.hidden_channels = [];
%     end
%     if isfield(N2N3_EEG,'psd'); N2N3_EEG = rmfield(N2N3_EEG, 'psd');end
%     for n = 1:N2N3_EEG.trials
%         [N2N3_EEG.psd.data(:,:,:,n),N2N3_EEG.psd.Hzbins] = psddata(squeeze(N2N3_EEG.data(:,:,n)),N2N3_EEG.srate,6,40,1); % only use fast for massive data
%       
%     end
%     close all;
%     plot_bands_spectra_EEG_PAMMS_Eisai(N2N3_EEG);
%     corr_ans = 0;
%     while ~corr_ans
%         finished = input('Finished Y/N? [Y]','s');
%         if ~ismember(finished, {'Y','N'})
%             corr_ans = 0;
%         else
%             corr_ans = 1;
%         end
%     end
% end
% 
% %save file data set
% tic
% if finished == 'Y'
    finishedfilename = [strrep(n2n3filename,'_n2n3_noarousals_fftstep_artifcorr.set','_n2n3_noarousals_fftstep_artifcorr_final'),'.set'];
    %           finishedfilename = [strrep(EEG.filename,'_n2n3_noarousals_artifcorr_fftstep.set','_n2n3_noarousals_artifcorr_fftstep_final'),'.set'];
    N2N3_EEG.datfile = finishedfilename;
    N2N3_EEG = pop_saveset(N2N3_EEG,'filename', finishedfilename, 'filepath',  N2N3_EEG.save_dir, 'check', 'on');
% end
disp(['... elapsed time ',num2str(toc/60),' saving set file']);

sprintf('run next section')
%% interpolate the removed channels
% ````````````````````````````````
%tic
%N2N3_EEG = eeg_interp(N2N3_EEG, N2N3_EEG.original_chanlocs);
%disp(['... elapsed time ',num2str(toc/60),' interpolating channels']);

%close all;

%if isfield(N2N3_EEG,'psdinterp'); N2N3_EEG = rmfield(N2N3_EEG, 'psdinterp');end
%for n = 1:N2N3_EEG.trials
    %[N2N3_EEG.psdinterp.data(:,:,:,n),N2N3_EEG.psdinterp.Hzbins] = psddata(squeeze(N2N3_EEG.data(:,:,n)),N2N3_EEG.srate,6,40,1); % only use fast for massive data
    %             EEG.trialstagenumber(n) = nanmean(EEG.psg.data(1,:,n));
    %             ei = find([EEG.event.epoch] == n);
    %             EEG.trialevents{n} = unique({EEG.event(ei).type});
%end

%N2N3_EEG.nbchan = EEG.nbchan; Uncomment if removing a channel 
%plot_bands_spectra_EEG_PAMMS_Eisai(N2N3_EEG);


%finalfilename = strrep(finishedfilename,'.set','_interp_alicecorrected.set');
%     finalfilename = strrep(EEG.filename,'.set','_interp.set');
%N2N3_EEG.setname  = [N2N3_EEG.setname,'interp'];
%N2N3_EEG = pop_saveset(N2N3_EEG, 'filename', finalfilename, 'filepath', N2N3_EEG.save_dir, 'check', 'on');

%figure('position',[680 575 630 525]);
%topoplot(cell2mat(N2N3_EEG.bad_channels),N2N3_EEG.original_chanlocs,'electrodes','labels',...
    %'shading','interp','style','blank','whitebk','on');
%cd(N2N3_EEG.save_dir);

%set(gcf,'color','w','paperpositionmode','auto');
%print(gcf,'-dpng','-r500',[ N2N3_EEG.save_dir,filesep,N2N3_EEG.filename,'.png']);

%plot_bands_spectra_EEG_PAMMS_Eisai(N2N3_EEG);
%sprintf('run next section')

 %% segment out N2 and N3 out as separate files
 cd (N2N3_EEG.save_dir)
load('COV207_N2N3_stagenames_afterfft.mat');
% cd(['D:\UWMpipeline\Scripts\temp_Loop\' subj_id '\' ])
N2_EEG = N2N3_EEG;
%cd 'C:' 
N2cleaned_filename = strrep(finishedfilename,'.set','_N2only.set');
N2only_cleaned = (strcmp(N2N3_stages_afterfft_snipped,'N2'));
swa_selectStagesEEGLAB_BAR(N2_EEG, N2only_cleaned,'N2only_cleaned',N2cleaned_filename);
N2only_cleaned = pop_saveset(N2N3_EEG, 'filename', N2cleaned_filename, 'check', 'on');
 

N3_EEG = N2N3_EEG;
N3cleaned_filename = strrep(finishedfilename,'.set','_N3only.set');
%cd 'C:'
N3only_cleaned = (strcmp(N2N3_stages_afterfft_snipped,'N3'));
swa_selectStagesEEGLAB_BAR(N3_EEG, N3only_cleaned,'N3only_cleaned',N3cleaned_filename);
N3only_cleaned = pop_saveset(N2N3_EEG, 'filename', N3cleaned_filename, 'check', 'on');

%end

% cd( N2N3_EEG.save_dir);
% %cd(['D:\UWMpipeline\Scripts\temp_Loop\' subj_id '\' ])
% NREMpsd = N2N3_EEG.psd;
% save NREMpsd NREMpsd;