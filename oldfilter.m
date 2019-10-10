%% Import the data from 'alldata.m'
% addpath (genpath ('C:\Users\taham\OneDrive - UHN\DBS_PD'))
% addpath(genpath('C:\Users\taham\Dropbox'))

% %this one for idir's PC
% addpath(genpath('C:\Users\User\Desktop\Taha\EEG_data'))

%% De-comment this for using alldata.mat
% % %load alldata;
% % inp=5; %input('please enter the number of trial. 1= s05, 2=s06, 3=STNon, 4=STNoff, 5=Paul 6=500msec:   ');
% % NumCol = alldata(inp).nbchan;
% % dpts=length(alldata(inp).data(:,1)); %this returns the total number of datapoints
% % 
% % data_Matrix = [alldata(inp).times alldata(inp).data];

%% 10-20 nsystem EEG (These lines remove headers from text data, if implemented)
fid = fopen('chenlab_file1_data.txt');
data = textscan(fid, '%*s %f %*[^\n]','HeaderLines',1);
fid = fclose(fid);
vec = data{1,1};
L = length(vec);
data_Matrix = zeros(L,NumCol);

      % These Lines make an array out of the tab delimited data
for i = 1:NumCol
fid = fopen('chenlab_file1_data.txt');
data = textscan(fid, [repmat('%*s',1,i-1), '%f', '%*[^\n]'],'HeaderLines',1);
fid = fclose(fid); 
data_Matrix(:,i) = data{1,1};
end

%% Plottings
figure; hold on,
plot(data_Matrix(:,1),data_Matrix(:,17)); %13->17

%% Find the peaks to select the individual trials
% --- Note: we should finally use the DBS onset for this purpose
dt = 0.05;% msec
Fs = 1/dt;% msec
sig = data_Matrix(30e3:40e3,17); % selected timeframe (set for prep 1, can be automatized later, e.g. with input)
sig = (sig - mean(sig));
figure; plot(sig)

a_max = max(sig);
Amp_th = a_max/10;
indx = find([sig;0]>=Amp_th); %[0;sig]<Amp_th & 

%% event finder (for prep1,2,5): Uses the event marking on the data
% dbsonset = alldata(inp).events.latency;
% en=1;
% for i=1:dbsonset(end)
%     epoch(en)=data_Matrix(dbsonset,:);
%     en=en+1;
% end

%% continue
indx_pks = zeros(length(indx),1);

for i = 2:length(indx)-1
    [pks,locs] = max(sig(indx(i)-5:indx(i)+5));
    indx_pks(i) = indx(i) - 10 - 1 + locs;
end
L_sel = 100; %add meaningful comment here NTS
sig_N = zeros(length(indx),L_sel+1);

for i = 2:length(indx)-10 %length(indx)
 %  sig_N(i,:) = sig(indx_pks(i)-L_sel/2:indx_pks(i)+L_sel/2)';
    sig_N(i,:) = sig(indx_pks(i)-5 : indx_pks(i)+L_sel-5)';
end

figure; plot(sig_N(1:end,:)','k')
hold on,
plot(mean(sig_N,1),'r','LineWidth',3)
title('Individual trials')
% %% Interpolation --> it should be based on the onset of DBS recorded by device
% tx = 0:length(sig_N)-1;
% tx_In = 0:0.05:length(sig_N)-1;
% sig_In = zeros(length(tx_In),76);
% for i = 1:76 %length(indx)
% %     sig_N(i,:) = sig(indx_pks(i)-L_sel/2:indx_pks(i)+L_sel/2)';
%    sig_In(:,i) = interp1(tx,sig_N(i,:)',tx_In);
% end
% figure; plot(sig_In,'k')
% hold on,
% plot(mean(sig_In,2),'r','LineWidth',3)


%% Crop the individual trials from their min, and normalize them
sig_test = zeros(size(sig_N));
d = zeros(length(sig_N),1);
for i = 1:50%size(sig_N,1)
    [pks,ind_min] = min(sig_N(i,:));
    test = -sig_N(i,ind_min:end);
    s_test = (test-mean(test))/max((test-mean(test)));
    d(1:length(s_test),i) = s_test';

end
figure; plot(d)
title('individual trials trimmed from their min')


%% Fitting B-spline to each individual pulse (purpose: to further remove the common signal (information) from each indiviual trial)
L_seg = 300; % samples --> this is for the artifact interval only
Des = mean(d,2);%
Ni = Des(4:L_seg);% ; zeros(length(Des)-200,1)];% signal to be smoothed (estimated by spline method)
numBin = 80;
tt = (1:length(Ni))*dt;
b = tt(1):numBin*dt:tt(end);% knots
        myI = eye(length(b));
        X = spline(b,myI,tt)*dt;
        ki = (X*X')^-1 * X*Ni;
        mu_i_ = X'*ki;

mu_i = [mu_i_;zeros(length(Des)-L_seg,1)]; % The smoothed signal   

figure; hold on,
plot(Des,'k')
plot(mu_i,'r')
legend('Average of all trials','smoothed estimate')


%% TBA:Cut 10 datapoints off the start of each trial



%% Adaptive Filtering
mu = 1;                % NLMS step size
offset = 10;           % NLMS offset
err_nlms = zeros(size(d)-3);
for i = 1:size(d,2)
    lms = dsp.LMSFilter('Length',20,'Method','Normalized LMS','StepSize',mu,'LeakageFactor',1);
    [y,err] = lms(d(4:end,i),mu_i);
    err_nlms(:,i) = err;
end

figure; plot(err_nlms(:,1),'k')
hold on,
plot(d(:,1),'r')
legend('Reconstructed signal(one trial)','The original trial')
figure; plot(err_nlms,'k')
title('Reconstructed signal')