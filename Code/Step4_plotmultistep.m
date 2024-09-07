file1 = 'out_sample_predictions';
file2 = 'combinedData.csv';

prediction = readtable(file1);
raw = readtable(file2);

actual_RK_out_sample = raw{end-365:end, 2};  
actual_RV_out_sample = raw{end-365:end, 3};  
actual_TSRV_out_sample = raw{end-365:end, 4};  
actual_RVpa_out_sample = raw{end-365:end, 5};  
predictions_RV_300 = prediction{:, 1:5};
predictions_RK = prediction{:, 6:10};
predictions_TSRV = prediction{:, 11:15};
predictions_RVpa = prediction{:, 16:20};


% step1
step1_RV = predictions_RV_300(:,1);
step1_RK = predictions_RK(:,1);
step1_TSRV = predictions_TSRV(:,1);
step1_RVpa = predictions_RVpa(:,1);


figure;
plot(actual_RV_out_sample, 'black', 'LineWidth', 1, 'DisplayName', 'Actual RV');
hold on;  
plot(step1_RV, 'r', 'LineWidth', 1, 'DisplayName', 'Predicted RV (Step 1)');
plot(step1_TSRV, 'yellow', 'LineWidth', 1, 'DisplayName', 'Predicted TSRV (Step 1)');
plot(step1_RK, 'green', 'LineWidth', 1, 'DisplayName', 'Predicted RK (Step 1)');
plot(step1_RVpa, 'blue', 'LineWidth', 1, 'DisplayName', 'Predicted RVpa (Step 1)');
title('out-of-sample Actual vs Predicted  for one-step-ahead');
xlabel('Time');
ylabel('Realized Volatility');
legend;
hold off;

saveas(gcf, 'com_multistep1.png');

% step2
step2_RV = predictions_RV_300(:,2);
step2_RK = predictions_RK(:,2);
step2_TSRV = predictions_TSRV(:,2);
step2_RVpa = predictions_RVpa(:,2);

figure;
plot(actual_RV_out_sample, 'black', 'LineWidth', 1, 'DisplayName', 'Actual RV');
hold on;  
plot(step2_RV, 'r', 'LineWidth', 1, 'DisplayName', 'Predicted RV (Step 2)');
plot(step2_TSRV, 'yellow', 'LineWidth', 1, 'DisplayName', 'Predicted TSRV (Step 2)');
plot(step2_RK, 'green', 'LineWidth', 1, 'DisplayName', 'Predicted RK (Step 2)');
plot(step2_RVpa, 'blue', 'LineWidth', 1, 'DisplayName', 'Predicted RVpa (Step 2)');
title('out-of-sample Actual vs Predicted  for two-step-ahead');
xlabel('Time');
ylabel('Realized Volatility');
legend;
hold off;

saveas(gcf, 'com_multistep2.png');

% step3
step3_RV = predictions_RV_300(:,3);
step3_RK = predictions_RK(:,3);
step3_TSRV = predictions_TSRV(:,3);
step3_RVpa = predictions_RVpa(:,3);

figure;
plot(actual_RV_out_sample, 'black', 'LineWidth', 1, 'DisplayName', 'Actual RV');
hold on;  
plot(step3_RV, 'r', 'LineWidth', 1, 'DisplayName', 'Predicted RV (Step 3)');
plot(step3_TSRV, 'yellow', 'LineWidth', 1, 'DisplayName', 'Predicted TSRV (Step 3)');
plot(step3_RK, 'green', 'LineWidth', 1, 'DisplayName', 'Predicted RK (Step 3)');
plot(step3_RVpa, 'blue', 'LineWidth', 1, 'DisplayName', 'Predicted RVpa (Step 3)');
title('out-of-sample Actual vs Predicted  for three-step-ahead');
xlabel('Time');
ylabel('Realized Volatility');
legend;
hold off;

saveas(gcf, 'com_multistep3.png');

% step4
step4_RV = predictions_RV_300(:,4);
step4_RK = predictions_RK(:,4);
step4_TSRV = predictions_TSRV(:,4);
step4_RVpa = predictions_RVpa(:,4);

figure;
plot(actual_RV_out_sample, 'black', 'LineWidth', 1, 'DisplayName', 'Actual RV');
hold on;  
plot(step4_RV, 'r', 'LineWidth', 1, 'DisplayName', 'Predicted RV (Step 4)');
plot(step4_TSRV, 'yellow', 'LineWidth', 1, 'DisplayName', 'Predicted TSRV (Step 4)');
plot(step4_RK, 'green', 'LineWidth', 1, 'DisplayName', 'Predicted RK (Step 4)');
plot(step4_RVpa, 'blue', 'LineWidth', 1, 'DisplayName', 'Predicted RVpa (Step 4)');
title('out-of-sample Actual vs Predicted  for four-step-ahead');
xlabel('Time');
ylabel('Realized Volatility');
legend;
hold off;

saveas(gcf, 'com_multistep4.png');

% step5
step5_RV = predictions_RV_300(:,5);
step5_RK = predictions_RK(:,5);
step5_TSRV = predictions_TSRV(:,5);
step5_RVpa = predictions_RVpa(:,5);

figure;
plot(actual_RV_out_sample, 'black', 'LineWidth', 1, 'DisplayName', 'Actual RV');
hold on;  
plot(step5_RV, 'r', 'LineWidth', 1, 'DisplayName', 'Predicted RV (Step 5)');
plot(step5_TSRV, 'yellow', 'LineWidth', 1, 'DisplayName', 'Predicted TSRV (Step 5)');
plot(step5_RK, 'green', 'LineWidth', 1, 'DisplayName', 'Predicted RK (Step 5)');
plot(step5_RVpa, 'blue', 'LineWidth', 1, 'DisplayName', 'Predicted RVpa (Step 5)');
title('out-of-sample Actual vs Predicted  for five-step-ahead');
xlabel('Time');
ylabel('Realized Volatility');
legend;
hold off;

saveas(gcf, 'com_multistep5.png');








