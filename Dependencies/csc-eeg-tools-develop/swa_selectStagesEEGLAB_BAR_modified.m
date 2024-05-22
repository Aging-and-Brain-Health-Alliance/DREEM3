% script to select stages from EEG1lab data
% ´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
% input
% load('C50_ds.set', '-mat')
% keepSamples = EEG1.swa_scoring.stages == 2;
% keepSamples(EEG1.swa_scoring.arousals) = false;


function swa_selectStagesEEG1LAB_BAR(EEG1, samples, setname,saveFile)
% check that samples is a logical array
if ~islogical(samples)
    fprintf(1, 'Error: second input should be a logical array');
end

% create the memory map to the data (not sure what this is doing) - It accelerates file data I/O by mapping the input data onto MATLAB's address space (like a mat array) //AD
tmp = memmapfile([EEG1.save_dir, '/', EEG1.datfile], 'Format', {'single', [EEG1.nbchan EEG1.pnts EEG1.trials], 'EEG1Data'}); %nbchan = no. of good channels; pnts = no. of samples/datapoints //AD
% the above formats the tmp data into a matrix called EEG1Data, which comprises (here) a 240 x 5547000 x 1 array into which the data in the specified filepath/filename location is stored
% select the desired samples for the memory mapped data
tmp1 = tmp.Data;
tmp2 = tmp1.EEG1Data;
newData = tmp2(:, samples); %in the vectorized EEG1Data, this picks out the samples labelled N2 and N3 but not their interseciton with arousal. 
% NOTE - this is also where the non-N2-N3-no-arousal is rejected BUT valid data is simply concatenated //AD

% adjust the EEG1 data
EEG1.times(~samples)     = []; %EEG1.times is just the second clock count of recording time. It has the same number of values as the no. of samples. Here, anything that isn't a valid N2/N3 sample has its timestamp DELETED - facilitating concatenation.
EEG1.pnts                = size(newData, 2); % recalculate the no. of datapoints based on whether they're N2/N3.
EEG1.data = newData;

% adjust the dataset name
EEG1.dataset = saveFile; 
EEG1.dataset(end-3: end)    = []; % cool way to get rid of a file extension and replace it //AD
EEG1.dataset = [EEG1.dataset, '.fdt'];
EEG1.datfile = EEG1.dataset;
EEG1.setname = setname;
%EEG1.filename = saveFile;
EEG1         = eeg_checkset(EEG1); % EEG1LAB function to test the consistency of data with the standard reference
% save the EEG1 struct
cd(EEG1.save_dir)
save(saveFile, 'EEG1', '-mat');%,'-v7.3');
display(['writing ',EEG1.save_dir saveFile]);

% write the new data file
fileID = fopen([EEG1.save_dir,EEG1.dataset],'w');
display(['writing ',EEG1.save_dir,EEG1.dataset]);
fwrite(fileID, newData,'single');
fclose(fileID);