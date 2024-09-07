% 
opts = detectImportOptions('Data_New.csv');
opts.VariableNames = {'Time', 'Price', 'Volume', 'Date'};
Data = readtimetable('Data_New.csv', opts);

unique_dates = unique(Data.Date);

% 
RVpa_table_BusinessTime = table('Size', [length(unique_dates), 2], 'VariableTypes', {'datetime', 'double'}, 'VariableNames', {'Date', 'RVpa_BusinessTime'});
RVpa_table_CalendarTime = table('Size', [length(unique_dates), 2], 'VariableTypes', {'datetime', 'double'}, 'VariableNames', {'Date', 'RVpa_CalendarTime'});

numDates = length(unique_dates);

% 
theta = 0.5; 

for date_idx = 1:numDates
    current_date = unique_dates(date_idx);
    day_data = Data(Data.Date == current_date, :);
    
    % 
    if height(day_data) < 2
        RVpa_table_BusinessTime{date_idx, 'Date'} = current_date;
        RVpa_table_BusinessTime{date_idx, 'RVpa_BusinessTime'} = NaN;
        RVpa_table_CalendarTime{date_idx, 'Date'} = current_date;
        RVpa_table_CalendarTime{date_idx, 'RVpa_CalendarTime'} = NaN;
        continue;
    end
    

    timeInSeconds = seconds(timeofday(day_data.Time)); 
    prices = day_data.Price; 
    
    valid_idx = ~isnan(prices);
    timeInSeconds = timeInSeconds(valid_idx);
    prices = prices(valid_idx);
    
    if length(prices) < 2
        RVpa_table_BusinessTime{date_idx, 'Date'} = current_date;
        RVpa_table_BusinessTime{date_idx, 'RVpa_BusinessTime'} = NaN;
        RVpa_table_CalendarTime{date_idx, 'Date'} = current_date;
        RVpa_table_CalendarTime{date_idx, 'RVpa_CalendarTime'} = NaN;
        continue;
    end
    
    % pre-average window(k_n)
    n = length(prices); 
    k_n = floor(theta * sqrt(n)); 
    
    % 
    timetype = 'seconds'; 
    samplinginterval = 1; 
    options = realized_options('Preaveraging'); 
    options.preAveragingKn = k_n; 
    
    % 
    try
        samplingtype_BusinessTime = 'BusinessTime';
        RVpa_BusinessTime = realized_preaveraged_variance(prices, timeInSeconds, timetype, samplingtype_BusinessTime, samplinginterval, options);
        RVpa_table_BusinessTime{date_idx, 'Date'} = current_date;
        RVpa_table_BusinessTime{date_idx, 'RVpa_BusinessTime'} = RVpa_BusinessTime; 
    catch ME
        disp(ME.message);
        RVpa_table_BusinessTime{date_idx, 'Date'} = current_date;
        RVpa_table_BusinessTime{date_idx, 'RVpa_BusinessTime'} = NaN;
    end
    
    % 
    try
        samplingtype_CalendarTime = 'CalendarTime';
        RVpa_CalendarTime = realized_preaveraged_variance(prices, timeInSeconds, timetype, samplingtype_CalendarTime, samplinginterval, options);
        RVpa_table_CalendarTime{date_idx, 'Date'} = current_date;
        RVpa_table_CalendarTime{date_idx, 'RVpa_CalendarTime'} = RVpa_CalendarTime; 
    catch ME
        disp(ME.message);
        RVpa_table_CalendarTime{date_idx, 'Date'} = current_date;
        RVpa_table_CalendarTime{date_idx, 'RVpa_CalendarTime'} = NaN;
    end
    
    % 
    fprintf('Progress: %.2f%% (%d/%d)\n', (date_idx/numDates)*100, date_idx, numDates);
end

figure;
subplot(2,1,1);
plot(RVpa_table_BusinessTime.Date, RVpa_table_BusinessTime.RVpa_BusinessTime);
xlabel('Date');
ylabel('RVpa BusinessTime');
title('RVpa using Business Time Sampling');
grid on;

subplot(2,1,2);
plot(RVpa_table_CalendarTime.Date, RVpa_table_CalendarTime.RVpa_CalendarTime);
xlabel('Date');
ylabel('RVpa CalendarTime');
title('RVpa using Calendar Time Sampling');
grid on;
saveas(gcf, 'RVpa_BusinessTime_CalendarTime.png');
writetable(RVpa_table_CalendarTime, 'RVpa.csv');