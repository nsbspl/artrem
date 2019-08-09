%% add paths
addpath (genpath ('C:\Users\taham\OneDrive - UHN\DBS_PD'))
addpath(genpath('C:\Users\taham\Dropbox'))

%% Import the data from central database in specific accessible variables for this code
load alldata % load the organized data recorded from patients and uploaded to the UHN OneDrive Cloud
trlnum=input('please enter the desired trial number. Cf. the patient listing on the DBS_PD folder');
NumCol = alldata(trlnum).nbchan; % 
