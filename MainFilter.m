%% Import the data

% path including file to be imported ** PLEASE SPECIFY
addpath("C:\Users\NeurotechUofT\Downloads\");
addpath("C:\Users\NeurotechUofT\Downloads\fieldtrip");
addpath(genpath("C:\Users\NeurotechUofT\Downloads\pwd"));
filename = 'Stop-signal_eDBS_Dec17_2018_s05-Block1.cnt';

%% importing file using FieldTrip
cfg                     = [];
cfg.dataset             = filename;
cfg.channel             = {'L2'};                                          % Please enter the desired channels here
cfg.continuous = 'yes';
raw = ft_preprocessing(cfg);


%% Import metadata

dpts = raw.sampleinfo(2); %total no. of datapoints
sampling_rate = raw.fsample; %2e4
dt = 1/sampling_rate; %0.05 msec

% helper variables:
NumCol = raw.sampleinfo(2);
inq_col = 1;
Nor = 1; %Normalize = 1 means that you will obtain normalized data for this program.
         % If you want the original scale, put Nor=0.

%displaying the imported metadata:
fprintf('this dataset is %s seconds long, equivalent to %d datapoints.\n',...
    num2str(length(raw.time{1,1})), dpts)
fprintf('the sampling rate is %s Hz.\n', num2str(sampling_rate))

%%
tbnds = [36, 36.5];
tdur = tbnds(2) - tbnds(1);
tsamp = tbnds.*sampling_rate;
sdur = tdur * sampling_rate;

time_secs = tbnds(1):dt:tbnds(2);%, sdur+1); %raw.time{1,1};
time_pts = tsamp(1):tsamp(2);

datamatrix = raw.trial{1,1}(tsamp(1):tsamp(2));
datamatrix = -(datamatrix - mean(datamatrix)); %makes the mean = 0 as well

figure; hold on,
plot(time_secs, datamatrix)
title('the imported section, mean = 0')
xlabel('time_{seconds}')
ylabel('v_{μvolts}')

%% Spike finder:

amp_max = max(datamatrix);
amp_th = amp_max/3;

[pks,locs] = findpeaks(datamatrix, 'MinPeakHeight', amp_th);

figure; hold on, findpeaks(datamatrix, 'MinPeakHeight', amp_th);
title('the imported section of the dataset, mean = 0')
xlabel('time_{seconds}')
ylabel('v_{μvolts}')

%% Cut the epochs:
% lstim = length(pk_idx);
% temp_len = 40; % to be set
% shift_left = 8;
% lstim = length(pks);

% GUIDE: cutter(self,locs,shift_left,lstim,temp_len)
eps = cutter(datamatrix,locs,8,length(pks),40);

%% Crop the individual trials from their min, and normalize them
sig_test = zeros(size(sig_N));
d = zeros(length(sig_N),1);

for i = 1:size(sig_N,1)
    [pks,ind_min] = min(sig_N(i,:));
    test = -sig_N(i,ind_min:end);
    if Nor 
        s_test = (test-mean(test))/max((test-mean(test)));
    else
        s_test = (test-mean(test));
    end
    d(1:length(s_test),i) = s_test';

end
figure; plot(time_vec_sel,d)
title('individual trials trimmed from their min')
% For each individual trial, this returns the the signal segment starting
% from the min. Then it is presented as the negative, so the starting point
% has the largest amplitude.

%% Fitting spline to each individual (purpose: to further remove the common signal (information) from each indiviual trial)
L_seg_ms = 10;
L_seg = L_seg_ms/dt; % samples --> this is for the artifact interval only
Des = mean(d,2);%
Ni = [Des(1:L_seg)];% ; zeros(length(Des)-200,1)];% signal to be smoothed (estimated by spline method)
numBin = 7;
tt = (1:length(Ni))*dt;
b = [tt(1):numBin*dt:tt(end)];% knots
        myI = eye(length(b));
        X = spline(b,myI,tt)*dt;
        ki = (X*X')^-1 * X*Ni;
        mu_i_ = X'*ki;

mu_i = [mu_i_;zeros(length(Des)-L_seg,1)]; % The smoothed signal   

figure; hold on,
plot(time_vec_sel,Des,'k')
plot(time_vec_sel,mu_i,'r')
legend('Average of all trials','smoothed estimate')

%% Adaptive Filtering
mu = 1;                % NLMS step size
offset = 10;           % NLMS offset
err_nlms = zeros(size(d));
for i = 1:size(d,2)
    lms = dsp.LMSFilter('Length',20,'Method','Normalized LMS','StepSize',mu,'LeakageFactor',1);
    [y,err] = lms(d(:,i),mu_i);
    err_nlms(:,i) = err;
end

figure; plot(time_vec_sel,err_nlms(:,1),'k')
hold on,
plot(time_vec_sel,d(:,1),'r')
legend('Reconstructed signal(one trial)','The original trial')
figure; plot(time_vec_sel,err_nlms,'k')
title('Reconstructed signal')
