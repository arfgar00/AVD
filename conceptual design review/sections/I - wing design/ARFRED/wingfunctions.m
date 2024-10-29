clc
clear
rho = -((12000-11887)/(12000-11500))*(0.336-0.311) + 0.336; % At 11887m AGL, linear regression between 11500 and 12000m altitude
u = 0.85 * 295.5; % Speed at given cruise
q = rho * u^2 / 2; % Dynamic pressure at given cruise
W_S = 8.144110697450234e+03; % W0/S
C_L_target = 0.7815; % Target C_L for angle of attack calculation

% Define the folder where your spreadsheet files are located
folder = '/Users/arfredgarcia/Documents/AVD/AVD GitHub Code/AVD/conceptual design review/functions/I - wing design/Airfoils';

% Get list of all CSV files in the folder
filePattern = fullfile(folder, '*.csv');  % For CSV files
csvFiles = dir(filePattern);

% Define the point where you want to find the gradient
x_point = 0.7815;

% Define the x range for interpolation
x_min = 0.7; % Lower limit for C_L fitting
x_max = 0.9; % Upper limit for C_L fitting

% Preallocate a structure to store results
results = struct('variableName', {}, 'LD', {}, 'alpha', {});

% Loop through each CSV file
for k = 1:length(csvFiles)
    % Get the full file name and path
    baseFileName = csvFiles(k).name;
    fullFileName = fullfile(folder, baseFileName);
    
    % Dynamically create variable name based on the file name (removing extension)
    [~, variableName, ~] = fileparts(baseFileName); % Extract file name without extension
    variableName = matlab.lang.makeValidName(variableName); % Ensure it's a valid variable name
    
    % Read the CSV file data (assume no headers)
    data = readtable(fullFileName);
    
    % Extract columns (1st column for alpha, 2nd column for C_L, and 3rd column for C_D)
    alpha = data{:, 1}; % Assuming the first column is angle of attack (alpha)
    C_L = data{:, 2};   % Assuming the second column is the lift coefficient (C_L)
    C_D = data{:, 3};   % Assuming the third column is the drag coefficient (C_D)

    % Step 1: Filter the data to only consider points between x_min and x_max
    valid_indices = (C_L >= x_min) & (C_L <= x_max); % Logical index for valid range
    C_L_filtered = C_L(valid_indices);
    C_D_filtered = C_D(valid_indices);
    alpha_filtered = alpha(valid_indices);
    
    % Step 2: Fit a cubic smoothing spline to the filtered data (C_L vs. C_D)
    smoothing_param = 0.9999;  % Adjust for smoother or closer fit (between 0 and 1)
    spline_fit_CD = csaps(C_L_filtered, C_D_filtered, smoothing_param);  % Fit the smoothing spline
    
    % Step 3: Evaluate the spline at 1000 points for a smooth curve
    C_L_fit = linspace(x_min, x_max, 1000); % Generate 1000 points within the range
    C_D_fit = fnval(spline_fit_CD, C_L_fit);  % Evaluate the spline at these points
    
    % Step 4: Create a single figure for both plots
    figure; % Create a new figure
    % Plot C_L vs. C_D on the first subplot
    subplot(2, 1, 1); % 2 rows, 1 column, first subplot
    plot(C_L, C_D, 'ro', 'DisplayName', 'All Data');  % Original data points
    hold on;
    plot(C_L_filtered, C_D_filtered, 'bo', 'DisplayName', 'Filtered Data');  % Filtered data points
    plot(C_L_fit, C_D_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Fitted Spline'); % Fitted spline
    
    % Step 5: Calculate the y value at x_point for the spline
    if x_point >= x_min && x_point <= x_max % Check if x_point is within the valid range
        C_D_at_x_point = fnval(spline_fit_CD, x_point); % Get the corresponding C_D value
        plot(x_point, C_D_at_x_point, 'g*', 'MarkerSize', 8, 'DisplayName', sprintf('Point at C_L = %.4f', x_point)); % Plot the point
    else
        fprintf('C_L = %.4f is outside the range for %s. Gradient not calculated.\n', x_point, variableName);
        continue; % Skip to the next iteration if out of range
    end
    
    % Additional plot formatting for C_D
    xlabel('C_L (Lift Coefficient)');
    ylabel('C_D (Drag Coefficient)');
    xlim([0, x_max]);
    ylim([0, 0.05]);
    legend('show');
    title(sprintf('Spline fitting for %s (C_L in [%.2f, %.2f])', variableName, x_min, x_max));
    
    % Step 6: Differentiate the spline to get the gradient (1st derivative)
    spline_derivative_CD = fnder(spline_fit_CD);  % Differentiate the spline
    
    % Step 7: Calculate the L/D at the specified x-point (C_L = 0.7815)
    if x_point >= x_min && x_point <= x_max % Check if x_point is within the valid range
        L_D_at_x = x_point / fnval(spline_fit_CD, x_point);
        % Store the results in the structure
        results(end + 1) = struct('variableName', variableName, 'LD', L_D_at_x, 'alpha', NaN); % Placeholder for alpha
    end
    
    % Step 8: Fit a cubic spline to alpha vs C_L
    spline_fit_alpha = csaps(C_L_filtered, alpha_filtered, smoothing_param);  % Fit the smoothing spline for alpha vs C_L
    
    % Step 9: Evaluate the spline for the alpha fit
    C_L_fit_alpha = linspace(min(C_L_filtered), max(C_L_filtered), 1000); % Generate 1000 points within the range
    alpha_fit_C_L = fnval(spline_fit_alpha, C_L_fit_alpha);  % Evaluate the spline at these points
    
    % Step 10: Plot alpha vs C_L on the second subplot
    subplot(2, 1, 2); % 2 rows, 1 column, second subplot
    plot(alpha, C_L, 'ro', 'DisplayName', 'All Data');  % Original data points for alpha vs C_L
    hold on;
    plot(alpha_filtered, C_L_filtered, 'bo', 'DisplayName', 'Filtered Data');  % Filtered data points
    plot(alpha_fit_C_L, C_L_fit_alpha, 'r-', 'LineWidth', 2, 'DisplayName', 'Fitted Spline (C_L vs Alpha)'); % Fitted spline for C_L vs alpha
    
    % Step 11: Specify the point C_L = 0.7815
    alpha_at_CL_target = fnval(spline_fit_alpha, C_L_target); % Get the corresponding alpha value
    plot(alpha_at_CL_target, C_L_target, 'g*', 'MarkerSize', 8, 'DisplayName', sprintf('Point at C_L = %.4f', C_L_target)); % Plot the point
    
    % Update the results structure with the corresponding alpha value
    results(end).alpha = alpha_at_CL_target; % Store alpha in the results
    
    % Additional plot formatting for alpha vs C_L
    xlabel('Angle of Attack (alpha)');
    ylabel('C_L (Lift Coefficient)');
    xlim([-10, 10]);
    ylim([-2, 2]);
    legend('show', 'Location', 'SouthEast');
    title(sprintf('Spline fitting for %s (Alpha vs C_L)', variableName));
    
    hold off; % Release the hold on the plot
end

% Sort the results based on the L/D values
[~, sort_indices] = sort([results.LD]); % Sort by L/D values
sorted_results = results(sort_indices); % Rearrange results according to sorted indices

% Display the sorted results
fprintf('L/D values at C_L = %.4f in order of magnitude:\n', x_point);
for i = 1:length(sorted_results)
    fprintf('L/D for %s: %.4f, Alpha: %.4f\n', sorted_results(i).variableName, sorted_results(i).LD, sorted_results(i).alpha);
end
