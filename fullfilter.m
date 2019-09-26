%% add paths

% the first two, for taha's laptop
addpath (genpath ('C:\Users\taham\OneDrive - UHN\DBS_PD'))
%addpath(genpath('C:\Users\taham\Dropbox'))

%this one for idir's PC
%addpath(genpath('C:\Users\User\Desktop\Taha\EEG_data'))
%

%% Import the data from central database to data_matrix [accessible variable for the purpose of this code]
load alldata % load the organized data recorded from patients and uploaded to the UHN OneDrive Cloud
trlnum=input('please enter the desired trial number. Cf. the patient listing on the DBS_PD folder  ');
cnl=input ('please add the number of ');
NumCol = alldata(trlnum).nbchan; % 

data_Matrix = [alldata(trlnum).times alldata(trlnum).data(:,cnl)];