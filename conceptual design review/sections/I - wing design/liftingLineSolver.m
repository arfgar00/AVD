clc;
clear;

% Step 1: Define wing and flight parameters for a Boeing 777 
b = 64.8;            % Span of the wing in meters
N = 101;              % Number of spanwise sections (must be odd for Simpson's rule)
y = linspace(-b/2, b/2, N);  % Spanwise positions from -b/2 to b/2
root_chord = 13.14;    % Root chord length in meters
tip_chord = 2.37;      % Tip chord length in meters


% Step 2: Define varying chord distribution along the span
chord = [linspace(tip_chord, root_chord, N/2), root_chord, linspace(root_chord, tip_chord, N/2)];

alpha = deg2rad(12.81);   % Angle of attack in radians
V_inf = 244;          % Freestream velocity in m/s (typical cruise speed)
rho = 1.225;          % Air density in kg/m^3
epsilon = 1e-13;       % Convergence criteria
D = 1e-2;             % Damping factor
a0 = 6.105;             % Lift curve slope for Boeing 777 airfoil in per radian (approx)

% Step 3: Initial guess for elliptical circulation distribution
Gamma_old = zeros(1, N); % Initialize Gamma_old to zeros
for i = 1:N
    if abs(y(i)) <= b/2
        Gamma_old(i) = sqrt(1 - (2*y(i)/b)^2); % Assign elliptical lift distribution
    else
        Gamma_old(i) = 0; % Set to zero outside the span
    end
end

figure 

plot(y,Gamma_old)

% Step 4 to 8: Iterative convergence loop
converged = false;
iteration = 0; % Counter to keep track of iteration
while ~converged
    iteration = iteration + 1;
    
    % Step 5: Calculate induced angle of attack using Lanchester-Prandtl lifting-line theory
    alpha_induced = zeros(1, N);
    dy = y(2) - y(1); % Interval size for Simpson's rule
    for i = 1:N
        % Apply Simpson's rule to calculate the induced downwash
        integral_sum = 0;
        for j = 1:N
            if i ~= j
                % Calculate distance
                distance = y(j) - y(i);
                
                % Check for singularity
                if abs(distance) < 1e-20 % Handling singularity case
                    % Use average of neighboring values
                    if i > 1 && i < N
                        distance = (y(j) - y(i-1) + y(j) - y(i+1)) / 2;
                    elseif i == 1
                        distance = (y(j) - y(i+1));
                    else
                        distance = (y(j) - y(i-1));
                    end
                end
                
                % Simpson's rule weighting
                if j == 1 || j == N
                    factor = 1; % Weight for first and last points
                elseif mod(j, 2) == 0
                    factor = 2; % Weight for even points
                else
                    factor = 4; % Weight for odd points
                end
    
                % Calculate the gradient of circulation at point i
                if i > 1 && i < N
                    grad_Gamma_i = (Gamma_old(i+1) - Gamma_old(i-1)) / (2 * dy); % Central difference
                elseif i == 1
                    grad_Gamma_i = (Gamma_old(2) - Gamma_old(1)) / dy; % Forward difference
                else
                    grad_Gamma_i = (Gamma_old(N) - Gamma_old(N-1)) / dy; % Backward difference
                end
                
                % Accumulate weighted contribution using the gradient
                integral_sum = integral_sum + factor / distance * grad_Gamma_i;
            end
        end
        integral_sum = (dy / 3) * integral_sum; % Apply Simpson's rule scaling
        alpha_induced(i) = integral_sum / (4 * pi * V_inf);
    end    
    % Step 6: Calculate effective angle of attack
    alpha_eff = alpha - alpha_induced;

    % Step 7: Obtain sectional lift coefficient
    Cl = a0 * alpha_eff;

    % Step 8: Calculate new circulation using Kutta-Joukowski theorem
    Gamma_new = Cl .* chord .* V_inf / 2;

    % Step 9: Check convergence
    diff = abs(Gamma_old - Gamma_new) ./ (abs(Gamma_new) ); 
    max_diff = max(diff); % Maximum relative difference for convergence check
    fprintf('Iteration %d: Max relative difference = %.20f\n', iteration, max_diff);
    
    if max_diff < epsilon
        converged = true;
        disp('Process is complete! Converged circulation distribution found.');
        break;
    end

    % Step 10: Update circulation input
    Gamma_old = Gamma_old + D * (Gamma_new - Gamma_old);

    if iteration > 100000
        break
    end
end

% Step 11: Calculate total circulation using Simpson's rule
total_circulation = (dy / 3) * (Gamma_old(1) + 4 * sum(Gamma_old(2:2:end)) + 2 * sum(Gamma_old(3:2:end))+ Gamma_old(end));

% Step 12: Calculate reference area (S)
S = (b) * mean(chord); % Reference area based on average chord length

% Step 13: Calculate sectional lift coefficient C_L
C_L = (2 / (V_inf * S)) * total_circulation;

% Step 14: Plot the final circulation distribution across the span
figure;
plot(y, Gamma_old, 'LineWidth', 1.5); % Plot circulation
xlabel('Spanwise Position (y)');
ylabel('Circulation (\Gamma)');
title('Converged Circulation Distribution for a Boeing 777');
grid on;

% Step 16: Plot the wing shape
figure;
hold on;

% Calculate the x-coordinates for the leading and trailing edges
x_leading_edge = zeros(1, N);  % Leading edge x-coordinates
x_trailing_edge = zeros(1, N); % Trailing edge x-coordinates
y_c = linspace(-b/2, b/2, N);   % Spanwise positions

for i = 1:N
    % Calculate the leading and trailing edge positions
    x_leading_edge(i) = chord(i) / 2;  % Leading edge x-coordinate
    x_trailing_edge(i) = -chord(i) / 2; % Trailing edge x-coordinate
end

% Create the coordinates for the wing shape
x_coords = [x_trailing_edge, flip(x_leading_edge)];
y_coords = [y_c, flip(y_c)];

% Plot the wing shape
fill(y_coords, x_coords, 'b', 'FaceAlpha', 0.5);

% Customize the plot
ylabel('Chord Position (m)');
xlabel('Spanwise Position (y)');
title('Wing Shape of Boeing 777');
grid on;
axis equal;
hold off;


% Step 16: Plot the final circulation distribution across the span
figure;
plot(y, Gamma_old, 'LineWidth', 1.5); % Plot circulation
xlabel('Spanwise Position (y)');
ylabel('Circulation (\Gamma)');
title('Converged Circulation Distribution for a Boeing 777');
grid on;

% Step 17: Output the total lift coefficient C_L
fprintf('Total Lift Coefficient (C_L): %.6f\n', C_L);
