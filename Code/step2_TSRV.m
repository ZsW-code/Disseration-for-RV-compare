% Load and preprocess data
opts = detectImportOptions('Data_New.csv');
opts.VariableNames = {'Time', 'Price', 'Volume', 'Date'};
Data = readtimetable('Data_New.csv', opts);

%file1 = ('omega_IV.csv'); 
%Data1 = readtable(file1);

%omega_squared = Data1.omega_squared;
%IV = Data1.IV;

unique_dates = unique(Data.Date);
numDates = length(unique_dates);

% Initialize result table
TSRV_table = table('Size', [numDates, 2], 'VariableTypes', {'datetime', 'double'}, 'VariableNames', {'Date', 'TSRV'});

% Set sparse interval (ensure this is defined)
sparse_interval = 1200; 

for date_idx = 1:numDates
    current_date = unique_dates(date_idx);
    day_data = Data(Data.Date == current_date, :);
    
    if height(day_data) < 2
        TSRV_table{date_idx, 'Date'} = current_date;
        TSRV_table{date_idx, 'TSRV'} = NaN;
        continue;
    end
    
    % Time and prices processing
    timeInSeconds = seconds(timeofday(day_data.Time)); 
    prices = day_data.Price; 
    
    valid_idx = ~isnan(prices);
    timeInSeconds = timeInSeconds(valid_idx);
    prices = prices(valid_idx);
    
    if length(prices) < 2
        TSRV_table{date_idx, 'Date'} = current_date;
        TSRV_table{date_idx, 'TSRV'} = NaN;
        continue;
    end
    
     % Get corresponding omega_squared and IV values
    omega_sq = omega_squared(date_idx);
    iv = IV(date_idx);
    
    % Calculate n
    n = length(prices);
    
    % Calculate optimal frequency 
    % the optimal subsample should be (sqrt(12)*omega^2/IV*n^2/3)

    subsamples = floor(sqrt(length(prices))); %Limited by computer performance

    % Calculate TSRV
    timetype = 'seconds'; 
    samplingtype = 'CalendarTime'; 
    samplinginterval = 1; 
    options = realized_options('Twoscale'); 
    
    try
        [rvts, rvtsss, rvtsd, rvtsssd, diagnostics] = realized_twoscale_variance(prices, timeInSeconds, timetype, samplingtype, samplinginterval, subsamples, options);
        TSRV_table{date_idx, 'Date'} = current_date;
        TSRV_table{date_idx, 'TSRV'} = rvtsd; 
    catch ME
        disp(ME.message);
        TSRV_table{date_idx, 'Date'} = current_date;
        TSRV_table{date_idx, 'TSRV'} = NaN;
    end
    
    fprintf('Progress: %.2f%% (%d/%d)\n', (date_idx / numDates) * 100, date_idx, numDates);
end

% Plot and save results
figure;
plot(TSRV_table.Date, TSRV_table.TSRV);
xlabel('Date');
ylabel('TSRV');
title('TSRV');
grid on;
saveas(gcf, 'TSRV.png');
writetable(TSRV_table, 'TSRV.csv');