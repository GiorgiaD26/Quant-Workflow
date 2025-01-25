clear all
close all
clc

% Load data from Excel
file_name = 'MSCI_Matlab.xlsx';
sheet_name = 'Sheet1';
data = readtable(file_name, 'Sheet', sheet_name);

% Extract dates
dates = data{:, 1};

% Convert dates to datetime format
if isnumeric(dates)
    dates = datetime(dates, 'ConvertFrom', 'excel');
else
    dates = datetime(dates, 'InputFormat', 'dd/MM/yyyy');
end

% Filter dates to get only those from April 30, 1974 to December 31, 2023
start_date = datetime(1974, 4, 30);
end_date = datetime(2023, 12, 31);
date_filter = (dates >= start_date) & (dates <= end_date);

% Apply filter to dates and prices
filtered_dates = dates(date_filter);
index_prices = data{date_filter, 2:end};

% Find columns without NaNs
nan_cols = any(isnan(index_prices), 1);
valid_indices = ~nan_cols;

% Select only indices with valid data (no NaNs)
filtered_labels = data.Properties.VariableNames(2:end); % Exclude the first column which contains dates
filtered_labels = filtered_labels(valid_indices); % Select only valid indices
filtered_prices = index_prices(:, valid_indices);

% Calculate monthly returns only for indices with valid data
filtered_returns = diff(log(filtered_prices)); % Calculate logarithmic returns (if data is prices)

% Calculate correlation matrix of returns
correlation_matrix = corr(filtered_returns, 'Rows', 'complete');

% Find correlations below 0.5 (considering symmetry)
correlation_matrix = triu(correlation_matrix, 1); % Consider only upper triangular part without diagonal
[row, col] = find(correlation_matrix < 0.5);

% Prepare data for table as columns
index1 = filtered_labels(row)';
index2 = filtered_labels(col)';
correlations = correlation_matrix(sub2ind(size(correlation_matrix), row, col));

% Create a mask to remove duplicates and correlations < 0.5
mask = row < col;
index1 = index1(mask);
index2 = index2(mask);
correlations = correlations(mask);

% Debugging: Check array dimensions
disp('Dimensions of index1:');
disp(size(index1));
disp('Dimensions of index2:');
disp(size(index2));
disp('Dimensions of correlations:');
disp(size(correlations));

% Check if array dimensions are equal
if length(index1) == length(index2) && length(index2) == length(correlations)
    % Create the table
    correlation_table = table(index1, index2, correlations, 'VariableNames', {'Index1', 'Index2', 'Correlation'});

    % Export table to a CSV file
    writetable(correlation_table, 'correlations_below_0.5.csv');

    % Display the table
    disp('Pairs of indices with correlation below 0.5:');
    disp(correlation_table);
else
    disp('Error: Array dimensions do not match.');
end

% Define sub-periods
periods = {
    datetime(1974, 4, 30), datetime(1984, 12, 31);
    datetime(1985, 1, 1), datetime(1995, 12, 31);
    datetime(1996, 1, 1), datetime(2006, 12, 31);
    datetime(2007, 1, 1), datetime(2017, 12, 31);
    datetime(2018, 1, 1), datetime(2023, 12, 31)
};

% Initialize a cell array for correlation matrices, labels, and sheet names
correlation_data = cell(size(periods, 1), 1);
sheet_names = cell(size(periods, 1), 1);

% Loop to calculate correlations for each sub-period
for i = 1:size(periods, 1)
    start_date = periods{i, 1};
    end_date = periods{i, 2};
    
    % Filter dates for the current period
    date_filter = (dates >= start_date) & (dates <= end_date);
    
    % Apply filter to dates and prices
    filtered_dates_period = dates(date_filter);
    filtered_prices_period = filtered_prices(date_filter, :);

    % Calculate monthly returns only for indices with valid data
    filtered_returns_period = diff(log(filtered_prices_period)); % Calculate logarithmic returns (if data is prices)

    % Calculate correlation matrix of returns
    correlation_matrix_period = corr(filtered_returns_period, 'Rows', 'complete');

    % Store correlation matrix and corresponding indices
    correlation_data{i}.Period = sprintf('%s_%s', datestr(start_date, 'yyyy-mm-dd'), datestr(end_date, 'yyyy-mm-dd'));
    correlation_data{i}.Indices = filtered_labels;
    correlation_data{i}.CorrelationMatrix = correlation_matrix_period;
    
    % Generate sheet name based on period
    sheet_names{i} = correlation_data{i}.Period;
end

% Create a new Excel file
excel_file = 'correlation_matrices_with_labels.xlsx';

% Write correlation matrices and index labels to separate Excel sheets
for i = 1:size(periods, 1)
    % Select current correlation matrix
    correlation_matrix = correlation_data{i}.CorrelationMatrix;
    indices_labels = correlation_data{i}.Indices;
    
    % Create a cell with index labels and correlation matrix
    % Insert index labels as row and column headers
    % Ensure indices_labels is a cell array of strings
    indices_labels = indices_labels';  % Transpose to match dimensions
    
    % Build data cell
    cell_data = cell(size(correlation_matrix, 1) + 1, size(correlation_matrix, 2) + 1);
    
    % Insert index labels as row and column headers
    cell_data(2:end, 1) = indices_labels;  % Left column labels
    cell_data(1, 2:end) = indices_labels;  % Top row labels
    
    % Insert correlation matrix into data cell
    cell_data(2:end, 2:end) = num2cell(correlation_matrix);
    
    % Write data cell to Excel sheet
    sheet_name = sheet_names{i};
    writecell(cell_data, excel_file, 'Sheet', sheet_name, 'Range', 'A1');
end

%%
% Calculate number of months in the dataset
num_months = size(filtered_returns, 1);

% Set lengths of moving windows in terms of months
window_3years = 36; % 3 years considering 12 months/year
window_5years = 60; % 5 years considering 12 months/year

% Preallocate matrices to store correlations
corr_3years = NaN(num_months - window_3years + 1, length(filtered_labels), length(filtered_labels));
corr_5years = NaN(num_months - window_5years + 1, length(filtered_labels), length(filtered_labels));

% Calculate correlations with moving windows
for i = 1:length(filtered_labels)
    for j = i+1:length(filtered_labels) % Use j = i+1:length(filtered_labels) to consider only unique pairs
        % Rolling window of 3 years
        for t = 1:num_months - window_3years + 1
            corr_3years(t, i, j) = corr(filtered_returns(t:t+window_3years-1, i), filtered_returns(t:t+window_3years-1, j));
        end
        
        % Rolling window of 5 years
        for t = 1:num_months - window_5years + 1
            corr_5years(t, i, j) = corr(filtered_returns(t:t+window_5years-1, i), filtered_returns(t:t+window_5years-1, j));
        end
    end
end

% Graph correlations over time for all unique asset pairs
num_assets = length(filtered_labels);
num_subplots = nchoosek(num_assets, 2); % Total number of subplots

% Set subplot grid dimensions
subplot_rows = 3; % Number of subplot rows per figure
subplot_cols = 3; % Number of subplot columns per figure
subplots_per_figure = subplot_rows * subplot_cols; % Number of subplots per figure

subplot_index = 1; % Current subplot index
figure_index = 1; % Current figure index

for i = 1:num_assets
    for j = i+1:num_assets % Use j = i+1:num_assets to consider only unique pairs
        % Get correlations for pair (i, j)
        corr_3years_plot = squeeze(corr_3years(:, i, j));
        corr_5years_plot = squeeze(corr_5years(:, i, j));
        
        % Calculate dates for plot in years
        plot_years_3years = year(filtered_dates(window_3years:num_months));
        plot_years_5years = year(filtered_dates(window_5years:num_months));
        
        % If necessary, create a new figure
        if mod(subplot_index - 1, subplots_per_figure) == 0
            if figure_index > 1
                % Save previous figure in JPG format
                saveas(gcf, sprintf('figure_%d.jpg', figure_index - 1));
                % Add legend to previous figure
                figure(figure_index - 1);
                legend('Rolling Window 3 years', 'Rolling Window 5 years', 'Location', 'northoutside');
            end
            figure;
            figure_index = figure_index + 1;
            subplot_counter = 1; % Reset subplot counter for new figure
        end
        
        % Draw subplot
        subplot(subplot_rows, subplot_cols, subplot_counter);
        plot(plot_years_3years, corr_3years_plot, 'b-', 'LineWidth', 1.5);
        hold on;
        plot(plot_years_5years, corr_5years_plot, 'r--', 'LineWidth', 1.5);
        title(sprintf('%s vs %s', filtered_labels{i}, filtered_labels{j}));
        xlabel('Years'); % Adjust x-axis label
        xlim([1974 2023])
        ylabel('Correlation');
        grid on;
        
        % Add trendline
        trendline_3years = polyfit(plot_years_3years, corr_3years_plot, 1);
        trendline_5years = polyfit(plot_years_5years, corr_5years_plot, 1);
        plot(plot_years_3years, polyval(trendline_3years, plot_years_3years), 'b-', 'LineWidth', 1);
        plot(plot_years_5years, polyval(trendline_5years, plot_years_5years), 'r--', 'LineWidth', 1);
        
        % Increment indices
        subplot_index = subplot_index + 1;
        subplot_counter = subplot_counter + 1;
    end
end

% Save the last figure in JPG format
saveas(gcf, sprintf('figure_%d.jpg', figure_index - 1));

% Add legend to the last figure
figure(figure_index - 1);
legend('Rolling Window 3 years', 'Rolling Window 5 years', 'Location', 'best');

% Additional settings for graph layout
sgtitle('Correlation Variation over Time for Unique Asset Pairs');
