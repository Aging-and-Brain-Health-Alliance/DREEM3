%% CONVERTS EDF TO SET FILE + SLEEP STAGES TO SAMPLES (for Eisai dataset)
close all
clear
clc
addpath(genpath('D:\UWMpipeline\Scripts\Preprocessing'))
addpath(genpath('D:\UWMpipeline\eeglab2021.0'))
%rmpath('D:\UWMpipeline\eeglab2021.0\plugins\Biosig3.7.5\biosig\maybe-missing\')
eeglab;

basePath = 'G:\Eisai\Clinilab';
EEG_dir = dir(basePath);

for subjectCnt = 33 : 37 %length(EEG_dir)
    subjectNum = EEG_dir(subjectCnt).name;
    sub_dir = dir([EEG_dir(subjectCnt).folder '\' subjectNum]);
    for sub_dir_cnt = 3 : 10 %length(sub_dir)
        
        for PSG = 2:7
            if strfind(sub_dir(sub_dir_cnt).name, ['PSG0' num2str(PSG)])
                destinationDir = dir([sub_dir(sub_dir_cnt).folder '\' sub_dir(sub_dir_cnt).name]);
                
                subj_id = ['Eisai_' subjectNum 'PSG0' num2str(PSG)];
                efd_File = dir(fullfile(destinationDir(1).folder, '*.edf'));
                sta_File = dir(fullfile(destinationDir(1).folder, '*.sta'));
                if ~isempty(efd_File) && ~isempty(sta_File)
                    temporary_save_directory=['D:\UWMpipeline\Scripts\temp_Loop_REM\' subj_id '\']
                    save_dir = temporary_save_directory;
                    mkdir(temporary_save_directory)
                    EEG_file = [destinationDir(1).folder '\' efd_File.name];
                    EEG = pop_biosig(EEG_file);
                    epoch_file = [destinationDir(1).folder '\' sta_File.name];
                    epoch=load(epoch_file);
                    epoch=epoch(:,2);%this column has the scores 7 being not scored, 0 being wake, 1&2&3&5 are sleep stages
                    epoch=epoch+1; %to be able to find the first 0 because find function searches for non-zero so from now on wake is 1, etc
                    Fs = EEG.srate;%frequency rate
                    first_epoch= find(epoch==1,1,'first');%first wake epch that now is 1 because we added 1
                    last_epoch= find (epoch==1,1,'last');%last wake epoch
                    sleep_stages=epoch(first_epoch:last_epoch);%in between first to last stage of scored sleep
                    
                    %% save and creat .set
                    filename = [subj_id '_Signals.set'];
                    
                    CHs={'eog','emg','f3','f4','c3','c4','o1','o2'};
                    CH_inds=[];
                    CHs_name={};
                    temp_chanlocs=EEG.chanlocs;
                    
                    cnt=1;
                    for ch1=1:length(temp_chanlocs)
                        for ch_ref=1:length(CHs)
                            if contains(lower(temp_chanlocs(ch1).labels),lower(CHs{ch_ref})) && ~contains(lower(temp_chanlocs(ch1).labels),'sp')
                                CHs_name{cnt}=temp_chanlocs(ch1).labels;
                                cnt=cnt+1;
                            end
                        end
                    end
                    
                    
                    EEG = pop_select( EEG, 'channel',CHs_name);
                    
                    
                    
                    EEG = pop_saveset(EEG, 'filename', filename, 'filepath', save_dir, 'check', 'on');
                    
                    %% clip the EDF based on sleep stages
                    first_data=(first_epoch-1)*Fs*30+1;
                    last_data=last_epoch*Fs*30;
                    
                    EEG = pop_select( EEG, 'point',[first_data last_data] );
                    EEG.sleepstage_lables=sleep_stages;
                    
                    %EEG.stagenames=
                    %to be added: an array 1x574000 showing the sleep stage at each sample EEG.something
                    %% mapping sleep stages to all samples of sleep
                    
                    %this loop is assinging a sleep stage to all the samples
                    sleep_stages2=repmat(sleep_stages,1,Fs*30);
                    sleep_stages2=reshape(sleep_stages2',size(sleep_stages2,1)*size(sleep_stages2,2),1)-1;
                    
                    
                    EEG.stagenames = sleep_stages2;
                    
                    theCells = cell(size(EEG.stagenames));
                    theCells(EEG.stagenames==0) = {'W'};
                    theCells(EEG.stagenames==1) = {'N1'};
                    theCells(EEG.stagenames==2) = {'N2'};
                    theCells(EEG.stagenames==3) = {'N3'};
                    theCells(EEG.stagenames==5) = {'REM'};
                    theCells(EEG.stagenames==6) = {'Arousal'};
                    theCells(EEG.stagenames==7) = {'NS'};
                    EEG.stagenames = theCells;
                    
                    %filtering to .35 to 35
                    EEG = pop_eegfiltnew(EEG, 'locutoff',59.5,'hicutoff',60.5,'revfilt',1);
                    EEG = pop_eegfiltnew(EEG, 'locutoff',0.35,'hicutoff',35);
                    
                    filename = sprintf('%s_Signals.set',subj_id);
                    EEG = pop_saveset(EEG, 'filename', filename, 'filepath', save_dir, 'check', 'on');
                end
            end
        end
    end
end

