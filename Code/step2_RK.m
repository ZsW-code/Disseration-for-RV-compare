%% Load data
opts = detectImportOptions('Data_New.csv');
opts.VariableNames = {'Time', 'Price', 'Volume', 'Date'};
Data = readtimetable('Data_New.csv', opts);

unique_dates = unique(Data.Date);

realizedKernelsTable = table('Size', [length(unique_dates), 2], 'VariableTypes', {'datetime', 'double'}, 'VariableNames', {'Date', 'RealizedKernel'});


numDates = length(unique_dates);
sparse_interval = 1200; 
omega_IV_Table = table('Size', [numDates, 3], 'VariableTypes', {'datetime', 'double', 'double'}, ...
                        'VariableNames', {'Date', 'omega_squared', 'IV'});

for date_idx = 1:numDates
    current_date = unique_dates(date_idx);
    day_data = Data(Data.Date == current_date, :);
    
    if height(day_data) < 2
        realizedKernelsTable{date_idx, 'Date'} = current_date;
        realizedKernelsTable{date_idx, 'RealizedKernel'} = NaN;
        continue;
    end
    
    prices = fillmissing(day_data.Price, 'previous');
    timeInSeconds = seconds(day_data.Time - day_data.Time(1));
    
    if length(prices) < 2
        realizedKernelsTable{date_idx, 'Date'} = current_date;
        realizedKernelsTable{date_idx, 'RealizedKernel'} = NaN;
        continue;
    end
    
    n = length(prices);
    q = max(1, floor(n / 720)); %make sure 2 minutes apart on average. 60*24/2=720
    log_returns = diff(log(prices));
    
    if any(isnan(log_returns)) || any(isinf(log_returns))
        disp(['Warning: NaN or Inf in log_returns on date ', datestr(current_date)]);
        realizedKernelsTable{date_idx, 'Date'} = current_date;
        realizedKernelsTable{date_idx, 'RealizedKernel'} = NaN;
        continue;
    end
    
    omega_squared_samples = zeros(q, 1);
    for i = 1:q
        subsample_indices = i:q:length(log_returns);
        subsampled_returns = log_returns(subsample_indices);
        RV_dense_i = sum(subsampled_returns.^2);
        n_i = length(subsampled_returns);
        omega_squared_samples(i) = RV_dense_i / (2 * n_i);
    end
    
    omega_squared = mean(omega_squared_samples);
    
    sparse_indices = 1:sparse_interval:length(log_returns);
    sparse_returns = log_returns(sparse_indices);
    IV = sum(sparse_returns.^2);
    
    min_IV_threshold = 1e-8;
    min_integrated_quarticity_threshold = 1e-12;
    
    if IV < min_IV_threshold
        IV = min_IV_threshold;
    end
    
    integrated_quarticity = IV^2;
    
    if integrated_quarticity < min_integrated_quarticity_threshold
        integrated_quarticity = min_integrated_quarticity_threshold;
    end
    omega_IV_Table{date_idx, 'Date'} = current_date;
    omega_IV_Table{date_idx, 'omega_squared'} = omega_squared;
    omega_IV_Table{date_idx, 'IV'} = integrated_quarticity;
   
    c_star = 3.5134;
    xi = sqrt(omega_squared / sqrt(integrated_quarticity));
    H_star = round(c_star * xi^(4/5) * n^(3/5)); 
     
    max_H_star = 5000; 
    if isnan(H_star) || H_star < 1
        H_star = 1; 
    elseif H_star > max_H_star
        H_star = max_H_star;
    end
    
    fprintf('Date: %s, H_star: %d, n: %d, q: %d, omega_squared: %.6f, IV: %.6f, integrated_quarticity: %.6f\n', ...
        datestr(current_date), H_star, n, q, omega_squared, IV, integrated_quarticity);
    
    try
        options = realized_options('Kernel'); 
        options.kernel = 'parzen';
        options.bandwidth = H_star; 
        options.useDebiasedNoise = true;
        
        [rk, rk_adjusted, diagnostics] = realized_kernel(prices, timeInSeconds, 'seconds', 'CalendarTime', 1, options);
        
        realizedKernelsTable{date_idx, 'Date'} = current_date;
        realizedKernelsTable{date_idx, 'RealizedKernel'} = rk;
        
    catch ME
        disp(['Error using MFE Toolbox on date ', datestr(current_date)]);
        disp(ME.message);
        realizedKernelsTable{date_idx, 'Date'} = current_date;
        realizedKernelsTable{date_idx, 'RealizedKernel'} = NaN;
    end
    
    fprintf('Progress: %.2f%% (%d/%d)\n', 100 * date_idx / numDates, date_idx, numDates);
end

figure;
plot(realizedKernelsTable.Date, realizedKernelsTable.RealizedKernel);
xlabel('Date');
ylabel('RK');
title('Realized Kernel over Time');
grid on;
saveas(gcf, 'RK.png'); 
writetable(realizedKernelsTable, 'RK.csv');

writetable(omega_IV_Table, 'omega_IV.csv');