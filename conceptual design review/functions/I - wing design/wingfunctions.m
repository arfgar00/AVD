clc
clear
rho = -((12000-11887)/(12000-11500))*(0.336-0.311) + 0.336; % at 11887m AGL, linear regression between 11500 and 12000m altitude
u = 0.85 * 295.5; % speed at given cruise
q = rho * u^2 / 2; % dynamic pressure at given cruise
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
x_min = 0.7;
x_max = 0.9;

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
    
    % Extract columns (1st column for alpha, 2nd column for x, and 3rd column for C_L)
    alpha = data{:, 1}; % Assuming the first column is angle of attack (alpha)
    x = data{:, 2}; % Assuming the second column is the independent variable (x)
    y = data{:, 3}; % Assuming the third column is the dependent variable (y), corresponding to C_L

    % Step 1: Filter the data to only consider points between x_min and x_max
    valid_indices = (x >= x_min) & (x <= x_max); % Logical index for valid range
    x_filtered = x(valid_indices);
    y_filtered = y(valid_indices);
    alpha_filtered = alpha(valid_indices); % Filtered alpha values

    % Step 2: Fit a cubic smoothing spline to the filtered data
    smoothing_param = 0.9999;  % Adjust for smoother or closer fit (between 0 and 1)
    spline_fit = csaps(x_filtered, y_filtered, smoothing_param);  % Fit the smoothing spline
    
    % Step 3: Evaluate the spline at 1000 points for a smooth curve
    x_fit = linspace(x_min, x_max, 1000); % Generate 1000 points within the range
    y_fit = fnval(spline_fit, x_fit);  % Evaluate the spline at these points
    
    % Step 4: Plot the original data and the fitted spline curve
    figure;
    plot(x, y, 'ro', 'DisplayName', 'All Data');  % Original data points within range
    hold on;
    plot(x_filtered, y_filtered, 'bo', 'DisplayName', 'Filtered Data');  % Filtered data points
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Fitted Spline'); % Fitted spline
    
    % Step 5: Calculate the y value at x_point for the spline
    if x_point >= x_min && x_point <= x_max % Check if x_point is within the valid range
        y_at_x_point = fnval(spline_fit, x_point); % Get the corresponding y value
        plot(x_point, y_at_x_point, 'g*', 'MarkerSize', 8, 'DisplayName', sprintf('Point at x = %.4f', x_point)); % Plot the point
    else
        fprintf('x = %.4f is outside the range for %s. Gradient not calculated.\n', x_point, variableName);
        continue; % Skip to the next iteration if out of range
    end
    
    % Additional plot formatting
    xlabel('C_L (x)');
    ylabel('C_d (y)');
    xlim([0,x_max]);
    ylim([0,0.05]);
    legend('show');
    title(sprintf('Spline fitting for %s (x in [0, 2.5])', variableName));
    hold off;
    
    % Step 6: Differentiate the spline to get the gradient (1st derivative)
    spline_derivative = fnder(spline_fit);  % Differentiate the spline
    
    % Step 7: Calculate the L/D at the specified x-point (x = 0.7815)
    if x_point >= x_min && x_point <= x_max % Check if x_point is within the valid range
        L_D_at_x = x_point / fnval(spline_fit, x_point);
        
        % Step 8: Fit the C_L vs. alpha data to find angle of attack for the target C_L
        alpha_spline_fit = csaps(alpha_filtered, y_filtered, smoothing_param); % Fit spline on alpha
        % Generate a smooth curve for angle of attack fitting
        alpha_fit = linspace(min(alpha_filtered), max(alpha_filtered), 1000);
        cl_fit = fnval(alpha_spline_fit, alpha_fit); % Evaluate fitted C_L for corresponding alpha

        % Step 9: Find the angle of attack corresponding to the target C_L
        if any(cl_fit >= C_L_target)
            alpha_at_CL_target = interp1(cl_fit, alpha_fit, C_L_target, 'linear', 'extrap'); % Linear interpolation
        else
            fprintf('No valid C_L values found for %s to calculate alpha at C_L = %.4f\n', variableName, C_L_target);
            alpha_at_CL_target = NaN; % Assign NaN if no valid C_L found
        end

        % Store the results in the structure
        results(end + 1) = struct('variableName', variableName, 'LD', L_D_at_x, 'alpha', alpha_at_CL_target);
    end
end

% Sort the results based on the L/D values
[~, sort_indices] = sort([results.LD]); % Sort by L/D values
sorted_results = results(sort_indices); % Rearrange results according to sorted indices

% Display the sorted results
fprintf('L/D values at x = %.4f in order of magnitude:\n', x_point);
for i = 1:length(sorted_results)
    fprintf('L/D for %s: %.4f, Alpha: %.4f degrees\n', sorted_results(i).variableName, sorted_results(i).LD, sorted_results(i).alpha);
end
