%% load data
file1 = 'RQ_5min.csv';
file2 = 'RV_VSP.csv';
file3 = 'RK.csv';
file4 = 'TSRV.csv';
file5 = 'RVpa.csv';

Data1 = readtable(file1);
Data2 = readtable(file2);
Data3 = readtable(file3);
Data4 = readtable(file4);
Data5 = readtable(file5);


Date = Data3{:, 1};       
RQ = Data1{Data1{:, 2} > 0, 2};   %filter no trading days   details is on step-1
RV = Data2{Data2{:, 300} > 0, 300}; %filter no trading days
RK = Data3{:, 2};         
TSRV = Data4{:, 2};      
RVpa = Data5{:, 2};      


combinedData = table(Date, RK, RV, TSRV, RVpa, RQ);

combinedData.Properties.VariableNames = {'Date', 'RealizedKernel', 'RV_300', 'TSRV', 'RVpa', 'RQ'};

writetable(combinedData, 'combinedData.csv');

%% drasmatic analysis
% Calculate mean, standard deviation, skewness, kurtosis, and log for each series
%RV
RV_mean = mean(RV);
RV_std = std(RV);
RV_skewness = skewness(RV);
RV_kurtosis = kurtosis(RV);
RV_min = min(RV);
RV_max = max(RV);
% log transform
RV_log = log(RV);
RV_log_mean = mean(RV_log);
RV_log_std = std(RV_log);
RV_log_skewness = skewness(RV_log);
RV_log_kurtosis = kurtosis(RV_log);
RV_log_min = min(RV_log);
RV_log_max = max(RV_log);

% Display RV statistics
disp('Statistical Analysis of RV:');
disp(['Mean: ', num2str(RV_mean)]);
disp(['Standard Deviation: ', num2str(RV_std)]);
disp(['Skewness: ', num2str(RV_skewness)]);
disp(['Kurtosis: ', num2str(RV_kurtosis)]);
disp(['Min: ', num2str(RV_min)]);
disp(['Max: ', num2str(RV_max)]);
disp(['Log Mean: ', num2str(RV_log_mean)]);
disp(['Log Standard Deviation: ', num2str(RV_log_std)]);
disp(['Log Skewness: ', num2str(RV_log_skewness)]);
disp(['Log Kurtosis: ', num2str(RV_log_kurtosis)]);
disp(['Log Min: ', num2str(RV_log_min)]);
disp(['Log Max: ', num2str(RV_log_max)]);
disp(' ');

%RK
RK_mean = mean(RK);
RK_std = std(RK);
RK_skewness = skewness(RK);
RK_kurtosis = kurtosis(RK);
RK_min = min(RK);
RK_max = max(RK);
% log transform
RK_log = log(RK);
RK_log_mean = mean(RK_log);
RK_log_std = std(RK_log);
RK_log_skewness = skewness(RK_log);
RK_log_kurtosis = kurtosis(RK_log);
RK_log_min = min(RK_log);
RK_log_max = max(RK_log);

% Display RK statistics
disp('Statistical Analysis of RK:');
disp(['Mean: ', num2str(RK_mean)]);
disp(['Standard Deviation: ', num2str(RK_std)]);
disp(['Skewness: ', num2str(RK_skewness)]);
disp(['Kurtosis: ', num2str(RK_kurtosis)]);
disp(['Min: ', num2str(RK_min)]);
disp(['Max: ', num2str(RK_max)]);
disp(['Log Mean: ', num2str(RK_log_mean)]);
disp(['Log Standard Deviation: ', num2str(RK_log_std)]);
disp(['Log Skewness: ', num2str(RK_log_skewness)]);
disp(['Log Kurtosis: ', num2str(RK_log_kurtosis)]);
disp(['Log Min: ', num2str(RK_log_min)]);
disp(['Log Max: ', num2str(RK_log_max)]);
disp(' ');



%TSRV
TSRV_mean = mean(TSRV);
TSRV_std = std(TSRV);
TSRV_skewness = skewness(TSRV);
TSRV_kurtosis = kurtosis(TSRV);
TSRV_min = min(TSRV);
TSRV_max = max(TSRV);
% log transform
TSRV_log = log(TSRV);
TSRV_log_mean = mean(TSRV_log);
TSRV_log_std = std(TSRV_log);
TSRV_log_skewness = skewness(TSRV_log);
TSRV_log_kurtosis = kurtosis(TSRV_log);
TSRV_log_min = min(TSRV_log);
TSRV_log_max = max(TSRV_log);

% Display TSRV statistics
disp('Statistical Analysis of TSRV:');
disp(['Mean: ', num2str(TSRV_mean)]);
disp(['Standard Deviation: ', num2str(TSRV_std)]);
disp(['Skewness: ', num2str(TSRV_skewness)]);
disp(['Kurtosis: ', num2str(TSRV_kurtosis)]);
disp(['Min: ', num2str(TSRV_min)]);
disp(['Max: ', num2str(TSRV_max)]);
disp(['Log Mean: ', num2str(TSRV_log_mean)]);
disp(['Log Standard Deviation: ', num2str(TSRV_log_std)]);
disp(['Log Skewness: ', num2str(TSRV_log_skewness)]);
disp(['Log Kurtosis: ', num2str(TSRV_log_kurtosis)]);
disp(['Log Min: ', num2str(TSRV_log_min)]);
disp(['Log Max: ', num2str(TSRV_log_max)]);
disp(' ');

%RVpa
RVpa_mean = mean(RVpa);
RVpa_std = std(RVpa);
RVpa_skewness = skewness(RVpa);
RVpa_kurtosis = kurtosis(RVpa);
RVpa_min = min(RVpa);
RVpa_max = max(RVpa);

% log transform
RVpa_log = log(RVpa);
RVpa_log_mean = mean(RVpa_log);
RVpa_log_std = std(RVpa_log);
RVpa_log_skewness = skewness(RVpa_log);
RVpa_log_kurtosis = kurtosis(RVpa_log);
RVpa_log_min = min(RVpa_log);
RVpa_log_max = max(RVpa_log);

% Display RVpa statistics
disp('Statistical Analysis of RVpa:');
disp(['Mean: ', num2str(RVpa_mean)]);
disp(['Standard Deviation: ', num2str(RVpa_std)]);
disp(['Skewness: ', num2str(RVpa_skewness)]);
disp(['Kurtosis: ', num2str(RVpa_kurtosis)]);
disp(['Min: ', num2str(RVpa_min)]);
disp(['Max: ', num2str(RVpa_max)]);
disp(['Log Mean: ', num2str(RVpa_log_mean)]);
disp(['Log Standard Deviation: ', num2str(RVpa_log_std)]);
disp(['Log Skewness: ', num2str(RVpa_log_skewness)]);
disp(['Log Kurtosis: ', num2str(RVpa_log_kurtosis)]);
disp(['Log Min: ', num2str(RVpa_log_min)]);
disp(['Log Max: ', num2str(RVpa_log_max)]);
disp(' ');

