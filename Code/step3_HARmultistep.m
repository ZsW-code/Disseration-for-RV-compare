file1 = 'combinedData.csv';
Data = readtable(file1);

RealizedKernel = log(Data{:, 2}); 
RV_300 = log(Data{:, 3});         
TSRV = log(Data{:, 4});           
RVpa = log(Data{:, 5});    
RQ = Data{:, 6};  
timeData = Data{:, 1};

n = length(RealizedKernel);

endDate = timeData(end);
startInSampleDate = endDate - years(4);  % in sample start
endInSampleDate = endDate - years(1);    % in sample end
startOutSampleDate = endDate - years(1); % out of sample end

% split sample
sampleInStartIndex = find(timeData >= startInSampleDate & timeData < startOutSampleDate, 1);  
sampleInEndIndex = find(timeData < startOutSampleDate, 1, 'last');  

sampleOutStartIndex = find(timeData >= startOutSampleDate, 1);  
sampleOutEndIndex = find(timeData <= endDate, 1, 'last');  

RealizedKernel_inSample = RealizedKernel(sampleInStartIndex:sampleInEndIndex);
RV_300_inSample = RV_300(sampleInStartIndex:sampleInEndIndex);
TSRV_inSample = TSRV(sampleInStartIndex:sampleInEndIndex);
RVpa_inSample = RVpa(sampleInStartIndex:sampleInEndIndex);


lag1_RK = lagmatrix(RealizedKernel_inSample, 1);
lag7_RK = mean(lagmatrix(RealizedKernel_inSample, 1:7), 2, 'omitnan');
lag30_RK = mean(lagmatrix(RealizedKernel_inSample, 1:30), 2, 'omitnan');

lag1_RV_300 = lagmatrix(RV_300_inSample, 1);
lag7_RV_300 = mean(lagmatrix(RV_300_inSample, 1:7), 2, 'omitnan');
lag30_RV_300 = mean(lagmatrix(RV_300_inSample, 1:30), 2, 'omitnan');

lag1_TSRV = lagmatrix(TSRV_inSample, 1);
lag7_TSRV = mean(lagmatrix(TSRV_inSample, 1:7), 2, 'omitnan');
lag30_TSRV = mean(lagmatrix(TSRV_inSample, 1:30), 2, 'omitnan');

lag1_RVpa = lagmatrix(RVpa_inSample, 1);
lag7_RVpa = mean(lagmatrix(RVpa_inSample, 1:7), 2, 'omitnan');
lag30_RVpa = mean(lagmatrix(RVpa_inSample, 1:30), 2, 'omitnan');


X_RK = [ones(length(lag1_RK), 1), lag1_RK, lag7_RK, lag30_RK];
X_RK = X_RK(30:end, :);  
y_RK = RealizedKernel_inSample(30:end);  

X_RV_300 = [ones(length(lag1_RV_300), 1), lag1_RV_300, lag7_RV_300, lag30_RV_300];
X_RV_300 = X_RV_300(30:end, :); 
y_RV_300 = RV_300_inSample(30:end);  

X_TSRV = [ones(length(lag1_TSRV), 1), lag1_TSRV, lag7_TSRV, lag30_TSRV];
X_TSRV = X_TSRV(30:end, :);  
y_TSRV = TSRV_inSample(30:end);  

X_RVpa = [ones(length(lag1_RVpa), 1), lag1_RVpa, lag7_RVpa, lag30_RVpa];
X_RVpa = X_RVpa(30:end, :);  
y_RVpa = RVpa_inSample(30:end);  

%% Estimation of HAR model using OLS
% RV_300 HAR with OLS
[b_RV_300, bint_RV_300, r_RV_300, rint_RV_300, stats_RV_300] = regress(y_RV_300, X_RV_300);

se_RV_300 = (bint_RV_300(:,2) - bint_RV_300(:,1)) / (2 * 1.96);  
t_values_RV_300 = b_RV_300 ./ se_RV_300;  
p_values_RV_300 = 2 * (1 - tcdf(abs(t_values_RV_300), length(y_RV_300) - length(b_RV_300)));  

fprintf('HAR Model Coefficients for RV_300:\n');
disp('Coefficient | Standard Error | t-value | p-value');
for i = 1:length(b_RV_300)
    fprintf('%12.4f | %14.4f | %8.4f | %8.4f\n', b_RV_300(i), se_RV_300(i), t_values_RV_300(i), p_values_RV_300(i));
end
fprintf('R-squared: %f\n', stats_RV_300(1));
fprintf('F-statistic: %f\n', stats_RV_300(2));
fprintf('F-test p-value: %f\n', stats_RV_300(3));
fprintf('Standard Error of Estimate: %f\n', stats_RV_300(4));


% RealizedKernel HAR with OLS
[b_RK, bint_RK, r_RK, rint_RK, stats_RK] = regress(y_RV_300, X_RK);

se_RK = (bint_RK(:,2) - bint_RK(:,1)) / (2 * 1.96);  
t_values_RK = b_RK ./ se_RK;  
p_values_RK = 2 * (1 - tcdf(abs(t_values_RK), length(y_RK) - length(b_RK)));  

fprintf('HAR Model Coefficients for RealizedKernel:\n');
disp('Coefficient | Standard Error | t-value | p-value');
for i = 1:length(b_RK)
    fprintf('%12.4f | %14.4f | %8.4f | %8.4f\n', b_RK(i), se_RK(i), t_values_RK(i), p_values_RK(i));
end
fprintf('R-squared: %f\n', stats_RK(1));
fprintf('F-statistic: %f\n', stats_RK(2));
fprintf('F-test p-value: %f\n', stats_RK(3));
fprintf('Standard Error of Estimate: %f\n', stats_RK(4));

% TSRV HAR with OLS
[b_TSRV, bint_TSRV, r_TSRV, rint_TSRV, stats_TSRV] = regress(y_RV_300, X_TSRV);

se_TSRV = (bint_TSRV(:,2) - bint_TSRV(:,1)) / (2 * 1.96);  
t_values_TSRV = b_TSRV ./ se_TSRV;  
p_values_TSRV = 2 * (1 - tcdf(abs(t_values_TSRV), length(y_TSRV) - length(b_TSRV)));  

fprintf('HAR Model Coefficients for TSRV:\n');
disp('Coefficient | Standard Error | t-value | p-value');
for i = 1:length(b_TSRV)
    fprintf('%12.4f | %14.4f | %8.4f | %8.4f\n', b_TSRV(i), se_TSRV(i), t_values_TSRV(i), p_values_TSRV(i));
end
fprintf('R-squared: %f\n', stats_TSRV(1));
fprintf('F-statistic: %f\n', stats_TSRV(2));
fprintf('F-test p-value: %f\n', stats_TSRV(3));
fprintf('Standard Error of Estimate: %f\n', stats_TSRV(4));

% RVpa HAR with OLS
[b_RVpa, bint_RVpa, r_RVpa, rint_RVpa, stats_RVpa] = regress(y_RV_300, X_RVpa);

se_RVpa = (bint_RVpa(:,2) - bint_RVpa(:,1)) / (2 * 1.96);  
t_values_RVpa = b_RVpa ./ se_RVpa;  
p_values_RVpa = 2 * (1 - tcdf(abs(t_values_RVpa), length(y_RVpa) - length(b_RVpa)));  

fprintf('HAR Model Coefficients for RVpa:\n');
disp('Coefficient | Standard Error | t-value | p-value');
for i = 1:length(b_RVpa)
    fprintf('%12.4f | %14.4f | %8.4f | %8.4f\n', b_RVpa(i), se_RVpa(i), t_values_RVpa(i), p_values_RVpa(i));
end
fprintf('R-squared: %f\n', stats_RVpa(1));
fprintf('F-statistic: %f\n', stats_RVpa(2));
fprintf('F-test p-value: %f\n', stats_RVpa(3));
fprintf('Standard Error of Estimate: %f\n', stats_RVpa(4));

y_pred_RV_300 = X_RV_300 * b_RV_300;
y_pred_RK = X_RK * b_RK;
y_pred_TSRV = X_TSRV * b_TSRV;
y_pred_RVpa = X_RVpa * b_RVpa;

% MSE in-sample
MSE_RV_300 = mean((exp(y_RV_300) - exp(y_pred_RV_300)).^2);
MSE_RK = mean((exp(y_RV_300) - exp(y_pred_RK)).^2);
MSE_TSRV = mean((exp(y_RV_300) - exp(y_pred_TSRV)).^2);
MSE_RVpa = mean((exp(y_RV_300) - exp(y_pred_RVpa)).^2);

RMSE_RV_300 = sqrt(MSE_RV_300);
RMSE_RK = sqrt(MSE_RK);
RMSE_TSRV = sqrt(MSE_TSRV);
RMSE_RVpa = sqrt(MSE_RVpa);


fprintf('MSE for RV_300: %f\n', MSE_RV_300);
fprintf('MSE for RealizedKernel: %f\n', MSE_RK);
fprintf('MSE for TSRV: %f\n', MSE_TSRV);
fprintf('MSE for RVpa: %f\n', MSE_RVpa);

fprintf('RMSE for RV_300: %f\n', RMSE_RV_300);
fprintf('RMSE for RealizedKernel: %f\n', RMSE_RK);
fprintf('RMSE for TSRV: %f\n', RMSE_TSRV);
fprintf('RMSE for RVpa: %f\n', RMSE_RVpa);

% Calculate QLIKE for each model
QLIKE_RV_300 = mean(    (exp(y_RV_300) ./ exp(y_pred_RV_300)) - log(exp(y_RV_300) ./ exp(y_pred_RV_300)) - 1     );
QLIKE_RK = mean((exp(y_RV_300) ./ exp(y_pred_RK)) - log(exp(y_RV_300) ./ exp(y_pred_RK)) - 1);
QLIKE_TSRV = mean((exp(y_RV_300) ./ exp(y_pred_TSRV)) - log(exp(y_RV_300) ./ exp(y_pred_TSRV)) - 1);
QLIKE_RVpa = mean((exp(y_RV_300) ./ exp(y_pred_RVpa)) - log(exp(y_RV_300) ./ exp(y_pred_RVpa)) - 1);

% Display QLIKE for each model
fprintf('QLIKE for RV_300: %f\n', QLIKE_RV_300);
fprintf('QLIKE for RealizedKernel: %f\n', QLIKE_RK);
fprintf('QLIKE for TSRV: %f\n', QLIKE_TSRV);
fprintf('QLIKE for RVpa: %f\n', QLIKE_RVpa);



%% rolling window forecast
 window = 365;
 n_steps = 5;  % 设置要预测的步数，例如5步预测

% 初始化存储样本期外的多步预测结果的矩阵
predictions_RV_300 = zeros(sampleOutEndIndex - sampleOutStartIndex + 1, n_steps);
predictions_RK = zeros(sampleOutEndIndex - sampleOutStartIndex + 1, n_steps);
predictions_TSRV = zeros(sampleOutEndIndex - sampleOutStartIndex + 1, n_steps);
predictions_RVpa = zeros(sampleOutEndIndex - sampleOutStartIndex + 1, n_steps);

 for t = sampleOutStartIndex:sampleOutEndIndex
     currentWindow_RV_300 = RV_300(t-window:t-1);  
     currentWindow_RK = RealizedKernel(t-window:t-1);  
     currentWindow_TSRV = TSRV(t-window:t-1);  
     currentWindow_RVpa = RVpa(t-window:t-1);  

     %lag realized measures
     lag1_RV_300 = lagmatrix(currentWindow_RV_300, 1);  
     lag7_RV_300 = mean(lagmatrix(currentWindow_RV_300, 1:7), 2, 'omitnan');  
     lag30_RV_300 = mean(lagmatrix(currentWindow_RV_300, 1:30), 2, 'omitnan');  

     lag1_RK = lagmatrix(currentWindow_RK, 1);  
     lag7_RK = mean(lagmatrix(currentWindow_RK, 1:7), 2, 'omitnan'); 
     lag30_RK = mean(lagmatrix(currentWindow_RK, 1:30), 2, 'omitnan');  

     lag1_TSRV = lagmatrix(currentWindow_TSRV, 1); 
     lag7_TSRV = mean(lagmatrix(currentWindow_TSRV, 1:7), 2, 'omitnan');  
     lag30_TSRV = mean(lagmatrix(currentWindow_TSRV, 1:30), 2, 'omitnan');  

     lag1_RVpa = lagmatrix(currentWindow_RVpa, 1);  
     lag7_RVpa = mean(lagmatrix(currentWindow_RVpa, 1:7), 2, 'omitnan'); 
     lag30_RVpa = mean(lagmatrix(currentWindow_RVpa, 1:30), 2, 'omitnan');  

     X_RV_300 = [ones(length(lag1_RV_300), 1), lag1_RV_300, lag7_RV_300, lag30_RV_300];
     X_RK = [ones(length(lag1_RK), 1), lag1_RK, lag7_RK, lag30_RK];
     X_TSRV = [ones(length(lag1_TSRV), 1), lag1_TSRV, lag7_TSRV, lag30_TSRV];
     X_RVpa = [ones(length(lag1_RVpa), 1), lag1_RVpa, lag7_RVpa, lag30_RVpa];

     y_RV_300 = currentWindow_RV_300(2:end);  
     y_RK = currentWindow_RK(2:end);
     y_TSRV = currentWindow_TSRV(2:end);
     y_RVpa = currentWindow_RVpa(2:end);

     % delete the first NaN
     X_RV_300 = X_RV_300(2:end, :);
     X_RK = X_RK(2:end, :);
     X_TSRV = X_TSRV(2:end, :);
     X_RVpa = X_RVpa(2:end, :);

     [b_RV_300, ~, ~, ~, ~] = regress(y_RV_300, X_RV_300);
     [b_RK, ~, ~, ~, ~] = regress(y_RV_300, X_RK);
     [b_TSRV, ~, ~, ~, ~] = regress(y_RV_300, X_TSRV);
     [b_RVpa, ~, ~, ~, ~] = regress(y_RV_300, X_RVpa);

     nextLag1_RV_300 = currentWindow_RV_300(end);  
     nextLag7_RV_300 = mean(currentWindow_RV_300(end-6:end));  
     nextLag30_RV_300 = mean(currentWindow_RV_300(end-29:end));  

     nextLag1_RK = currentWindow_RK(end);  
     nextLag7_RK = mean(currentWindow_RK(end-6:end));  
     nextLag30_RK = mean(currentWindow_RK(end-29:end));  

     nextLag1_TSRV = currentWindow_TSRV(end);  
     nextLag7_TSRV = mean(currentWindow_TSRV(end-6:end));  
     nextLag30_TSRV = mean(currentWindow_TSRV(end-29:end));  
     
     nextLag1_RVpa = currentWindow_RVpa(end);  
     nextLag7_RVpa = mean(currentWindow_RVpa(end-6:end));  
     nextLag30_RVpa = mean(currentWindow_RVpa(end-29:end));  


     for step = 1:n_steps
        
        next_X_RV_300 = [1, nextLag1_RV_300, nextLag7_RV_300, nextLag30_RV_300];
        next_X_RK = [1, nextLag1_RK, nextLag7_RK, nextLag30_RK];
        next_X_TSRV = [1, nextLag1_TSRV, nextLag7_TSRV, nextLag30_TSRV];
        next_X_RVpa = [1, nextLag1_RVpa, nextLag7_RVpa, nextLag30_RVpa];
        
        predicted_RV_300 = next_X_RV_300 * b_RV_300;
        predicted_RK = next_X_RK * b_RK;
        predicted_TSRV = next_X_TSRV * b_TSRV;
        predicted_RVpa = next_X_RVpa * b_RVpa;
      

        predictions_RV_300(t - sampleOutStartIndex + 1, step) = predicted_RV_300;
        predictions_RK(t - sampleOutStartIndex + 1, step) = predicted_RK;
        predictions_TSRV(t - sampleOutStartIndex + 1, step) = predicted_TSRV;
        predictions_RVpa(t - sampleOutStartIndex + 1, step) = predicted_RVpa;
      

        nextLag1_RV_300 = predicted_RV_300; 
        nextlag7_RV_300 = ((lag7_RV_300 * 7) - currentWindow_RV_300(end-6) + predicted_RV_300) / 7;
        nextlag30_RV_300 = ((lag30_RV_300 * 30) - currentWindow_RV_300(end-29) + predicted_RV_300) / 30;
        
        nextLag1_RK = predicted_RK;
        nextlag7_RK = ((lag7_RK * 7) - currentWindow_RK(end-6) + predicted_RK) / 7;
        nextlag30_RK = ((lag30_RK * 30) - currentWindow_RK(end-29) + predicted_RK) / 30;
        
        nextLag1_TSRV = predicted_TSRV;
        nextlag7_TSRV = ((lag7_TSRV * 7) - currentWindow_TSRV(end-6) + predicted_TSRV) / 7;
        nextlag30_TSRV = ((lag30_TSRV * 30) - currentWindow_TSRV(end-29) + predicted_TSRV) / 30;
        
        nextLag1_RVpa = predicted_RVpa;   
        nextlag7_RVpa = ((lag7_RVpa * 7) - currentWindow_RVpa(end-6) + predicted_RVpa) / 7;
        nextlag30_RVpa = ((lag30_RVpa * 30) - currentWindow_RVpa(end-29) + predicted_RVpa) / 30;
    end

 end

 % QLIKE and MSE 

 qlike_RV_300 = zeros(sampleOutEndIndex - sampleOutStartIndex + 1, n_steps);
 qlike_RK = zeros(sampleOutEndIndex - sampleOutStartIndex + 1, n_steps);
 qlike_TSRV = zeros(sampleOutEndIndex - sampleOutStartIndex + 1, n_steps);
 qlike_RVpa = zeros(sampleOutEndIndex - sampleOutStartIndex + 1, n_steps);

 avg_qlike_RV_300 = zeros(1, n_steps);
 avg_qlike_RK = zeros(1, n_steps);
 avg_qlike_TSRV = zeros(1, n_steps);
 avg_qlike_RVpa = zeros(1, n_steps);

 mse_RV_300 = zeros(1, n_steps);
 mse_RK = zeros(1, n_steps);
 mse_TSRV = zeros(1, n_steps);
 mse_RVpa = zeros(1, n_steps);

 rmse_RV_300 = zeros(1, n_steps);
 rmse_RK = zeros(1, n_steps);
 rmse_TSRV = zeros(1, n_steps);
 rmse_RVpa = zeros(1, n_steps);

for t = sampleOutStartIndex:sampleOutEndIndex
    
    if t + n_steps - 1 <= sampleOutEndIndex  
        
        for step = 1:n_steps
            
           
            actual_RV = exp(RV_300(t + step - 1));  
            actual_RK = exp(RV_300(t + step - 1));  
            actual_TSRV = exp(RV_300(t + step - 1));  
            actual_RVpa = exp(RV_300(t + step - 1)); 
          
 
            prediction_RV = exp(predictions_RV_300(t - sampleOutStartIndex + 1, step));
            prediction_RK = exp(predictions_RK(t - sampleOutStartIndex + 1, step));
            prediction_TSRV = exp(predictions_TSRV(t - sampleOutStartIndex + 1, step));
            prediction_RVpa = exp(predictions_RVpa(t - sampleOutStartIndex + 1, step));
            
            
            % QLIKE
            qlike_RV_300(t - sampleOutStartIndex + 1, step) = (actual_RV / prediction_RV) - log(actual_RV / prediction_RV) - 1;
            avg_qlike_RV_300(step) = mean(qlike_RV_300(:, step), 'omitnan');
            
            qlike_RK(t - sampleOutStartIndex + 1, step) = (actual_RK / prediction_RK) - log(actual_RK / prediction_RK) - 1;
            avg_qlike_RK(step) = mean(qlike_RK(:, step), 'omitnan');
            
            qlike_TSRV(t - sampleOutStartIndex + 1, step) = (actual_TSRV / prediction_TSRV) - log(actual_TSRV / prediction_TSRV) - 1;
            avg_qlike_TSRV(step) = mean(qlike_TSRV(:, step), 'omitnan');
           
            qlike_RVpa(t - sampleOutStartIndex + 1, step) = (actual_RVpa / prediction_RVpa) - log(actual_RVpa / prediction_RVpa) - 1;
            avg_qlike_RVpa(step) = mean(qlike_RVpa(:, step), 'omitnan');

            mse_RV_300(step) = mean((actual_RV - prediction_RV).^2, 'omitnan');
            mse_RK(step) = mean((actual_RV - prediction_RK).^2, 'omitnan');
            mse_TSRV(step) = mean((actual_RV - prediction_TSRV).^2, 'omitnan');
            mse_RVpa(step) = mean((actual_RV - prediction_RVpa).^2, 'omitnan');
            
            rmse_RV_300(step) = sqrt(mse_RV_300(step));
            rmse_RK(step) = sqrt(mse_RK(step));
            rmse_TSRV(step) = sqrt(mse_TSRV(step));
            rmse_RVpa(step) = sqrt(mse_RVpa(step));

        end
    end
end

% out of sample RMSE
fprintf('Out-of-Sample RMSE Results:\n');
fprintf('Step | RMSE (RV_300) | RMSE (RK) | RMSE (TSRV) | RMSE (RVpa)\n');
for step = 1:n_steps
    fprintf('%4d | %14.6f | %12.6f | %14.6f | %12.6f\n', step, ...
            rmse_RV_300(step), rmse_RK(step), rmse_TSRV(step), rmse_RVpa(step));
end

% out of sample QLIKE
fprintf('\nOut-of-Sample QLIKE Results:\n');
fprintf('Step | QLIKE (RV_300) | QLIKE (RK) | QLIKE (TSRV) | QLIKE (RVpa)\n');
for step = 1:n_steps
    fprintf('%4d | %14.6f | %12.6f | %14.6f | %12.6f\n', step, ...
            avg_qlike_RV_300(step), avg_qlike_RK(step), avg_qlike_TSRV(step), avg_qlike_RVpa(step));
end


%% DM test
% 初始化存储DM检验结果的矩阵
dm_results_RK = zeros(n_steps, 2);  % DM statistics and p-values for RV_300 vs RK
dm_results_TSRV = zeros(n_steps, 2);  % DM statistics and p-values for RV_300 vs TSRV
dm_results_RVpa = zeros(n_steps, 2);  % DM statistics and p-values for RV_300 vs RVpa

% 定义损失函数为MSE
dfun = @(x, y) x.^2 - y.^2;

% 计算预测误差和DM检验
for step = 1:n_steps
    
    % 计算 RV_300 和其他测度的预测误差
    errors_RV_300 = exp(RV_300(sampleOutStartIndex:sampleOutEndIndex)) - exp(predictions_RV_300(:, step));
    errors_RK = exp(RV_300(sampleOutStartIndex:sampleOutEndIndex)) - exp(predictions_RK(:, step));
    errors_TSRV = exp(RV_300(sampleOutStartIndex:sampleOutEndIndex)) - exp(predictions_TSRV(:, step));
    errors_RVpa = exp(RV_300(sampleOutStartIndex:sampleOutEndIndex)) - exp(predictions_RVpa(:, step));

    % DM检验: RV_300 vs RK
    [dm_stat_RK, p_value_RK] = dmtest_modified(errors_RV_300, errors_RK,dfun, step);
    dm_results_RK(step, :) = [dm_stat_RK, p_value_RK];

    % DM检验: RV_300 vs TSRV
    [dm_stat_TSRV, p_value_TSRV] = dmtest_modified(errors_RV_300, errors_TSRV,dfun,step);
    dm_results_TSRV(step, :) = [dm_stat_TSRV, p_value_TSRV];

    % DM检验: RV_300 vs RVpa
    [dm_stat_RVpa, p_value_RVpa] = dmtest_modified(errors_RV_300, errors_RVpa,dfun,step);
    dm_results_RVpa(step, :) = [dm_stat_RVpa, p_value_RVpa];

end

% 展示结果
fprintf('DM Test Results (using MSE as loss function):\n');
fprintf('Step | DM Stat (RK) | p-value (RK) | DM Stat (TSRV) | p-value (TSRV) | DM Stat (RVpa) | p-value (RVpa)\n');
for step = 1:n_steps
    fprintf('%4d | %12.4f | %12.4f | %15.4f | %14.4f | %15.4f | %14.4f\n', step, ...
            dm_results_RK(step, 1), dm_results_RK(step, 2), ...
            dm_results_TSRV(step, 1), dm_results_TSRV(step, 2), ...
            dm_results_RVpa(step, 1), dm_results_RVpa(step, 2));
end



%% model confidence set
n_models = 4;
n_samples = sampleOutEndIndex - sampleOutStartIndex + 1;  
n_steps = 5;  

% MCS (all default)
alpha = 0.1;  
B = 1000;      
W = 12;        

for step = 1:n_steps
    %loss function
    losses = zeros(n_samples, n_models);
    losses(:, 1) = qlike_RV_300(:, step);
    losses(:, 2) = qlike_RK(:, step);
    losses(:, 3) = qlike_TSRV(:, step);
    losses(:, 4) = qlike_RVpa(:, step);
    
    
    [includedR, pvalsR, excludedR, includedSQ, pvalsSQ, excludedSQ] = mcs(losses, alpha, B, W, 'STATIONARY');
    
    
    fprintf('Step %d:\n', step);
    disp('Included models using R method:');
    disp(includedR);
    
    disp('P-values using R method:');
    disp(pvalsR);
    
    disp('Excluded models using R method:');
    disp(excludedR);
    
    disp('Included models using SQ method:');
    disp(includedSQ);
    
    disp('P-values using SQ method:');
    disp(pvalsSQ);
    
    disp('Excluded models using SQ method:');
    disp(excludedSQ);
    
    fprintf('------------------------------------\n');
end

%% save predictions
predictions_out_sample = [exp(predictions_RV_300), exp(predictions_RK), exp(predictions_TSRV), exp(predictions_RVpa)];

column_names = {'RV_300_Step1', 'RV_300_Step2', 'RV_300_Step3', 'RV_300_Step4', 'RV_300_Step5', ...
                'RK_Step1', 'RK_Step2', 'RK_Step3', 'RK_Step4', 'RK_Step5', ...
                'TSRV_Step1', 'TSRV_Step2', 'TSRV_Step3', 'TSRV_Step4', 'TSRV_Step5', ...
                'RVpa_Step1', 'RVpa_Step2', 'RVpa_Step3', 'RVpa_Step4', 'RVpa_Step5'};


predictions_table = array2table(predictions_out_sample, 'VariableNames', column_names);
writetable(predictions_table, 'out_sample_predictions.csv');



