%% Import the data from 'alldata.m'
%load alldata;
inp = 5; %input('please enter the number of trial. 1= s05, 2=s06, 3=STNon, 4=STNoff, 5=TMS-EEG:   ');
NumCol = alldata(inp).nbchan;
dpts = length(alldata(inp).data(:,1)); %this returns the total number of datapoints
inq_col = 1;
dt = 0.05;% msec
Fs = 1/dt;% msec
Nor = 1; %Normalize = 1 means that you will obtain normalized data for this program.
         % If you want the original scale, put Nor=0.
data_Matrix = [alldata(inp).times, alldata(inp).data(:,inq_col)];
% for inp=5, total length = 135s

%% 10-20 nsystem EEG (These lines remove headers from text data, if implemented)
% NumCol=2;
% fid = fopen('chenlab_file1_data.txt');
% data = textscan(fid, '%*s %f %*[^\n]','HeaderLines',1);
% fid = fclose(fid);
% vec = data{1,1};
% L = length(vec);
% data_Matrix = zeros(L,NumCol);
% 
%       % This makes an array out of the tab delimited data
% for i = 1:NumCol
% fid = fopen('chenlab_file1_data.txt');
% data = textscan(fid, [repmat('%*s',1,i-1), '%f', '%*[^\n]'],'HeaderLines',1);
% fid = fclose(fid); 
% data_Matrix(:,i) = data{1,1};
% end
% inq_col = 2;

%% Plottings
figure; hold on,
plot(data_Matrix(:,1),data_Matrix(:,2));


%% Find the peaks to select the individual trials
% --- Note: we should finally use the DBS onset for this purpose
sig = data_Matrix(35e3:250e3,2); % selected timeframe (set for prep 1, can be automatized later, e.g. with input)
% 35e3:250e3 is 10.75s, you can change this.
sig = -(sig - mean(sig));
figure; plot(data_Matrix(35e3:250e3,1),sig)

a_max = max(sig);
Amp_th = a_max/10;
indx = find([0;sig]<Amp_th & [sig;0]>=Amp_th);

%% event finder (for prep1,2,5): Uses the event marking on the data
% dbsonset = alldata(inp).events.latency;
% en=1;
% for i=1:dbsonset(end)
%     epoch(en)=data_Matrix(dbsonset,:);
%     en=en+1;
% end

%% continue
indx_pks = zeros(length(indx),1);

for i = 1:length(indx)
    a = 10; % 2a is the length of searching interval for the peaks.                  
    [pks,locs] = max(sig(indx(i)-a:indx(i)+a));
    indx_pks(i) = indx(i) - a - 1 + locs;
end

L_sel_ms = 0.2e3;% each individual selected segment length is 200ms, you can change this.
L_sel = L_sel_ms/dt; % number of time points in the above segment length.
time_vec_sel = linspace(0,L_sel_ms,L_sel+1);
sig_N = zeros(length(indx),L_sel+1);

for i = 1:length(indx)
    sig_N(i,:) = sig(indx_pks(i) - 20 : indx_pks(i)+L_sel - 20)'; 
    % selected individual segment with length L_sel_ms starts from "1ms to the left of the peak"
end

figure; plot(time_vec_sel, sig_N(1:end,:)','k')
hold on,
plot(time_vec_sel, mean(sig_N,1),'r','LineWidth',3) % mean of all individual trials at each time point;
                                      % the red line is thicker
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