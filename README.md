# DREEM3
This repository contains the functions necessary for processing the eeg data derived from DREEM3 headbands.
To perform analyses on your data derive3d from DREEM3 device, you can run the SM_main_eeg_preprocess.m file. 

#Prior to running this function, you need:
1. Have eeglab installed. (and have its path added to your MATLAB)
2. Have statistics toolbox installed 
The inputs for this function are: 
1. subjects id
2. source path: the path to your data (subject folder)
3. dep_path: the path you have saved the dependency functions 
(on default, it should be in ...\DREEM3\Dependencies)

# What does SM_main_eeg_preprocess.m do?

1. Converts the .edf file to .set file. which is a structure that can be read by eeglab toolbox.
2. Adds the lables to the sleep stages.
3. Saves the .set file with sleep stages 
4. Pops up a visual control for the first time(you can check the overal quality of the whole data-REM and NREM- here)
5. Creates a .set structure for N2N3 (which represent the NREM sleep)
6. Detects and delete the bad signals using "wispic" method as implemented in eeglab
7. visually check the data after artifact rejection only for NREM sleep
8. check if you have even number of clips in deleting the bad signals





