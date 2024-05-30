% script to select stages from eeglab data
% ´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
% input
% load('C50_ds.set', '-mat')
% keepSamples = EEG.swa_scoring.stages == 2;
% keepSamples(EEG.swa_scoring.arousals) = false;


function swa_selectStagesEEGLAB_BAR(EEG, samples, setname,saveFile)
% check that samples is a logical array
if ~islogical(samples)
    fprintf(1, 'Error: second input should be a logical array');
end

% create the memory map to the data (not sure what this is doing) - It accelerates file data I/O by mapping the input data onto MATLAB's address space (like a mat array) //AD
tmp = memmapfile(fullfile(EEG.save_dir, EEG.datfile), 'Format', {'single', [EEG.nbchan EEG.pnts EEG.trials], 'eegData'}); %nbchan = no. of good channels; pnts = no. of samples/datapoints //AD
% the above formats the tmp data into a matrix called eegData, which comprises (here) a 240 x 5547000 x 1 array into which the data in the specified filepath/filename location is stored
% select the desired samples for the memory mapped data
tmp1 = tmp.Data;
tmp2 = tmp1.eegData;
newData = tmp2(:, samples); %in the vectorized eegData, this picks out the samples labelled N2 and N3 but not their interseciton with arousal. 
% NOTE - this is also where the non-N2-N3-no-arousal is rejected BUT valid data is simply concatenated //AD

% adjust the EEG data
EEG.times(~samples)     = []; %EEG.times is just the second clock count of recording time. It has the same number of values as the no. of samples. Here, anything that isn't a valid N2/N3 sample has its timestamp DELETED - facilitating concatenation.
EEG.pnts                = size(newData, 2); % recalculate the no. of datapoints based on whether they're N2/N3.
EEG.data = newData;
% adjust the dataset name
EEG.dataset = saveFile; 
EEG.dataset(end-3: end)    = []; % cool way to get rid of a file extension and replace it //AD
EEG.dataset = [EEG.dataset, '.fdt'];
EEG.datfile = EEG.dataset;
EEG.setname = setname;
%EEG.filename = saveFile;
EEG         = eeg_checkset(EEG); % EEGLAB function to test the consistency of data with the standard reference
% save the EEG struct
save(fullfile(EEG.save_dir,saveFile), 'EEG', '-mat');%,'-v7.3');
display(['writing ',fullfile(EEG.save_dir,EEG.dataset)]);

% write the new data file
fileID = fopen([EEG.save_dir,saveFile],'w');
display(['writing ',fullfile(EEG.save_dir,EEG.dataset)]);
fwrite(fileID, EEG.data,'single');
fclose(fileID);