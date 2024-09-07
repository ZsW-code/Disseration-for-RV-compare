%% Load data and preparation
opts = detectImportOptions('bitstampUSD.csv');
opts.VariableNames = {'Time', 'Price', 'Volume'};

data_org = readtable('bitstampUSD.csv', opts);
data_org.Time = datetime(data_org.Time, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
data_org.Date = dateshift(data_org.Time, 'start', 'day');

% Test negative volume
negativeVolumes = data_org.Volume < 0;
negativeData = data_org(negativeVolumes, :);
data_org.Volume = abs(data_org.Volume);

%plot
figure; 
yyaxis left;
plot(data_org.Time, data_org.Price, 'b', 'LineWidth', 2); 
ylabel('Price'); 
xlabel('Time');
ylim([min(data_org.Price) max(data_org.Price)]);

yyaxis right;
plot(data_org.Time, data_org.Volume, 'red'); 
ylabel('Volume'); 
ylim([0 max(data_org.Volume)*1.1])
title('Price and Volume Over Time(Tick-by-Tick)');
saveas(gcf, 'price_volume_plot_from2012.png');
disp("load done")

%% Define function for trimmed statistics
function [trimmed_mean, trimmed_std] = trimmed_stats(window, trim_percent)
    
    sorted_window = sort(window);
    trim_num = round(trim_percent * length(sorted_window));
    trimmed_window = sorted_window(trim_num+1:end-trim_num);
    trimmed_mean = mean(trimmed_window);
    trimmed_std = std(trimmed_window);
end

%% Identify outliers using trimmed mean and standard deviation
trim_percent = 0.01; 
prices = data_org.Price;
is_outlier = false(size(prices)); 

unique_dates = unique(data_org.Date);
for date_idx = 1:length(unique_dates)
    current_date = unique_dates(date_idx);
    date_indices = data_org.Date == current_date;
    date_prices = prices(date_indices);
    
    
    [trimmed_mean, trimmed_std] = trimmed_stats(date_prices, trim_percent);
    
    
    for i = find(date_indices)'
        if abs(prices(i) - trimmed_mean) > 5 * trimmed_std
            is_outlier(i) = true;
        end
    end
end


data_org.Outlier = is_outlier;
outliers = data_org(data_org.Outlier, :);
numOutliers = height(outliers);
data_clean = data_org(~data_org.Outlier, :);

figure;
subplot(2, 1, 1);  
plot(data_org.Time, data_org.Price, 'b.');  
hold on;
scatter(outliers.Time, outliers.Price, 'ro', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'SizeData', 20);  
hold off;
title('Price Time Series with Outliers Highlighted(Tick-by Tick)');
xlabel('Time');
ylabel('Price');
legend('Normal Data', 'Outliers');

subplot(2, 1, 2);  
plot(data_clean.Time, data_clean.Price, 'b.');  
title('Price Time Series without Outliers(Tick-by Tick)');
xlabel('Time');
ylabel('Price');
legend('Clean Data');
saveas(gcf, 'outliers.png');


% Test duplicate timestamp
[uniqueTimes, ia, ic] = unique(data_clean.Time);
duplicateIndices = setdiff(1:length(data_clean.Time), ia); 

duplicatedTimes = data_clean.Time(duplicateIndices);

% Check
if ~isempty(duplicatedTimes)
    disp('Exist duplicate timestamps');
else
    disp('No duplicate timestamps');
end

% Aggregate data
aggregatedVolume = accumarray(ic, data_clean.Volume, [], @sum);
weightedPriceSum = accumarray(ic, data_clean.Price .* data_clean.Volume, [], @sum);
aggregatedPrice = weightedPriceSum ./ aggregatedVolume;

% Create new table
Data = timetable(uniqueTimes, aggregatedPrice, aggregatedVolume, 'VariableNames', {'Price', 'Volume'});
Data = Data(Data.Volume > 0, :);
disp("pre done")

%% Initial plotting
figure; 
yyaxis left;
plot(Data.uniqueTimes, Data.Price, 'b', 'LineWidth', 2); 
ylabel('Price'); 
xlabel('Time');
ylim([min(Data.Price) max(Data.Price)]);

yyaxis right;
plot(Data.uniqueTimes, Data.Volume, 'red'); 
ylabel('Volume'); 
ylim([0 max(Data.Volume)*1.1])
title('Price and Volume Over Time after aggregating(Tick-by-Tick)');
grid on;
hold on
saveas(gcf, 'Data_after_aggregate.png');


% rearrange the time period(number of transcations > 1000)
cutoffdate = datetime(2014,1,1,'TimeZone','UTC');
Data_New = Data(Data.Properties.RowTimes >= cutoffdate, :);
Data_New.Date = dateshift(Data_New.uniqueTimes, 'start', 'day');  
%numbers of outliers after 2014.1.1
numOutliers_New = height(outliers(outliers.Date >= cutoffdate, :));

%% Identify and record no trading days
startDate = min(Data_New.Date);
endDate = max(Data_New.Date);
allDates = dateshift(startDate, 'start', 'day'):days(1):dateshift(endDate, 'start', 'day');
existingDates = unique(Data_New.Date);
missingDates = setdiff(allDates, existingDates);
noTradingDays = table(missingDates, 'VariableNames', {'NoTradingDates'});


figure; 
yyaxis left;
plot(Data_New.uniqueTimes, Data_New.Price, 'b', 'LineWidth', 2); 
ylabel('Price'); 
xlabel('Time');
ylim([min(Data_New.Price) max(Data_New.Price)]);

yyaxis right;
plot(Data_New.uniqueTimes, Data_New.Volume, 'red'); 
ylabel('Volume'); 
ylim([0 max(Data_New.Volume)*1.1])
title('Price and Volume Over Time(filtered Tick-by-Tick)');
grid on;

writetimetable(Data_New, 'Data_New.csv');
saveas(gcf, 'price_volume_comparison.png');

%% VSP and 5-Minute RQ Calculation
Data_New.Date = dateshift(Data_New.uniqueTimes, 'start', 'day');
minDate = min(Data_New.Date);
Data_New.DayNumber = days(Data_New.Date - minDate) + 1;

Tmax = max(Data_New.DayNumber);
RV_tick = zeros(Tmax, 1);
RV_VSP = zeros(Tmax, 3600);
RQ_5min = zeros(Tmax, 1); % Initialize RQ array

for i = 1:Tmax
    dayData = Data_New(Data_New.DayNumber == i, :);
    if isempty(dayData)
        disp(['Day ', num2str(i), ' has no trading data. Skipping...']);
        continue;
    end
    
    % Log prices and times
    ptemp = log(dayData.Price);
    ttemp = seconds(timeofday(dayData.Properties.RowTimes));
    
    % Calculate tick-level realized variance
    RV_tick(i) = sum(diff(ptemp).^2);
    
    % Calculate 5-minute intervals
    interval = 300; % 300 seconds = 5 minutes
    tgrid = 0:interval:86400;
    
    % Use histcounts to bin the timestamps into the grid intervals
    [counts, edges] = histcounts(ttemp - 0.01, tgrid);
    pind = cumsum(counts);
    
    % Ensure pind indices are within bounds
    pind(pind == 0) = 1;
    pind(pind > length(ptemp)) = length(ptemp);
    
    % Align prices with the grid
    p_grid = ptemp(pind);
    p_grid = [ptemp(1); p_grid]; % Add the opening prices
    
    % Forward fill NaN values in p_grid
    nan_indices = isnan(p_grid);
    for k = find(nan_indices)'
        if k > 1
            p_grid(k) = p_grid(k-1);
        else
            p_grid(k) = p_grid(find(~nan_indices, 1));
        end
    end
    
    % Calculate log returns for 5-minute intervals
    log_returns_5min = diff(p_grid);
    
    % Calculate Realized Quarticity (RQ) for the day
    if any(isnan(log_returns_5min)) || isempty(log_returns_5min)
        RQ_5min(i) = NaN;
    else
        RQ_5min(i) = sum(log_returns_5min.^4);
    end
    
    % Calculate RV_VSP for various intervals
    intervals = 1:3600;
    for m = 1:3600
        interval = intervals(m);
        tgrid = 0:interval:86400;
        
        % Use histcounts to bin the timestamps into the grid intervals
        [counts, edges] = histcounts(ttemp - 0.01, tgrid);
        pind = cumsum(counts);
        
        % Ensure pind indices are within bounds
        pind(pind == 0) = 1;
        pind(pind > length(ptemp)) = length(ptemp);
        
        % Align prices with the grid
        p_grid = ptemp(pind);
        p_grid = [ptemp(1); p_grid];
        
        % Forward fill NaN values in p_grid
        nan_indices = isnan(p_grid);
        for k = find(nan_indices)'
            if k > 1
                p_grid(k) = p_grid(k-1);
            else
                p_grid(k) = p_grid(find(~nan_indices, 1));
            end
        end
        
        % Calculate log returns and realized variance
        log_returns = diff(p_grid);
        
        if any(isnan(log_returns)) || isempty(log_returns)
            RV_VSP(i, m) = NaN;
        else
            RV_VSP(i, m) = sum(log_returns.^2);
        end
    end
    disp(['Processed day: ', num2str(i)]);
end

% Save the RQ results to a CSV file
RQ_table = table((1:Tmax)', RQ_5min, 'VariableNames', {'DayNumber', 'RQ_5min'});
writetable(RQ_table, 'RQ_5min.csv');

writematrix(RV_VSP, 'RV_VSP.csv');

mean_VSP = mean(RV_VSP); 

figure;
plot(1:3600, mean_VSP, 'LineWidth', 2);
xlabel('Sampling Interval');
xlim([0 3600])
ylabel('Average Realized Volatility');
title('Volatility Signature Plot');
grid on;
xline(300, 'label', 'A Viable Sampling Interval');
xticks([0  300 500 1000 1500 2000 2500 3000 3500 3600]);
xticklabels({'0', '300', '500', '1000', '1500', '2000', '2500', '3000', '3500'});
saveas(gcf, 'VSP.png');

RV_300 = RV_VSP(:,300);
RV_300(RV_300 == 0) = [];

log_RV_300 = log(RV_300);
figure;
histogram(log_RV_300, 30);
xlabel('Log Realized Volatility');
ylabel('Frequency');
title('Distribution of Log Realized Volatility (300-second Intervals)');


%Mean: -6.9286
%Standard Deviation: 1.1426
%Skewness: 0.097369
%Kurtosis: 3.6293

RV_300_table = table(existingDates,RV_300,'VariableNames', {'Date', 'RV_300'});
figure;
plot(RV_300_table.Date, RV_300_table.RV_300);
xlabel('Date');
ylabel('RV_300');
title('Realized Volatility (5-min Intervals) over Time');
grid on;
saveas(gcf, 'RV_300.png');
