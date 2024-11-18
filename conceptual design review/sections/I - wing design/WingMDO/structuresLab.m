% Constants for strain gauge conversion
Vin = 15;                % Input voltage to the strain gauge bridge (V)
S = 2.09;                % Strain gauge factor
G = 1000;                % Amplifier gain
L = 500;                 % Lever arm length (mm) for torque calculation
load_cases_front = [10, 20, 30, 40, 50] * 0.453592 * 9.81;  % Load increments for front spar (N)
load_cases_rear = [12, 18, 27, 36, 48] * 0.453592 * 9.81;   % Load increments for rear spar (N)

% Example dataset for strain gauge readings (All 9 gauges)
strain_gauge_data_front = [
    0.10, 0.15, 0.20, 0.25, 0.30;  % Gauge C (X-axis)
    0.20, 0.30, 0.40, 0.50, 0.60;  % Gauge D (45°)
    0.15, 0.25, 0.35, 0.45, 0.55;  % Gauge E (Y-axis)
    0.30, 0.40, 0.50, 0.60, 0.70;  % Gauge F (X-axis)
    0.40, 0.50, 0.60, 0.70, 0.80;  % Gauge G (45°)
    0.25, 0.35, 0.45, 0.55, 0.65;  % Gauge H (Y-axis)
    0.35, 0.45, 0.55, 0.65, 0.75;  % Gauge I (X-axis)
    0.45, 0.55, 0.65, 0.75, 0.85;  % Gauge J (45°)
    0.50, 0.60, 0.70, 0.80, 0.90;  % Gauge K (Y-axis)
];

% Example dataset for strain gauge readings (All 9 gauges)
strain_gauge_data_rear = [
    0.10, 0.15, 0.20, 0.25, 0.30;  % Gauge C (X-axis)
    0.20, 0.30, 0.40, 0.50, 0.60;  % Gauge D (45°)
    0.15, 0.25, 0.35, 0.45, 0.55;  % Gauge E (Y-axis)
    0.30, 0.40, 0.50, 0.60, 0.70;  % Gauge F (X-axis)
    0.40, 0.50, 0.60, 0.70, 0.80;  % Gauge G (45°)
    0.25, 0.35, 0.45, 0.55, 0.65;  % Gauge H (Y-axis)
    0.35, 0.45, 0.55, 0.65, 0.75;  % Gauge I (X-axis)
    0.45, 0.55, 0.65, 0.75, 0.85;  % Gauge J (45°)
    0.50, 0.60, 0.70, 0.80, 0.90;  % Gauge K (Y-axis)
];

% Number of strain gauges and load cases
num_gauges_rear = size(strain_gauge_data_rear, 1);
num_cases_rear = size(strain_gauge_data_rear, 2);
num_gauges_front = size(strain_gauge_data_front, 1);
num_cases_front = size(strain_gauge_data_front, 2);

% Initialize an array for strains
strains_front = zeros(num_gauges_front, num_cases_front);
strains_rear = zeros(num_gauges_front, num_cases_front);

% Step 1: Convert strain gauge voltage to strain
for i = 1:num_cases_front
    for j = 1:num_gauges_front
        V_out = strain_gauge_data_front(j, i);  % Voltage reading for each load case
        strains_front(j, i) = (4 * V_out) / (Vin * S * G);  % Convert voltage to strain
    end
end

for i = 1:num_cases_rear
    for j = 1:num_gauges_rear
        V_out = strain_gauge_data_rear(j, i);  % Voltage reading for each load case
        strains_rear(j, i) = (4 * V_out) / (Vin * S * G);  % Convert voltage to strain
    end
end

% Step 2: Group strains into rosettes
strains_CDE_front = strains_front(1:3, :);  % Group CDE (C: X-axis, D: 45°, E: Y-axis)
strains_FGH_front = strains_front(4:6, :);  % Group FGH (F: X-axis, G: 45°, H: Y-axis)
strains_IJK_front = strains_front(7:9, :);  % Group IJK (I: X-axis, J: 45°, K: Y-axis)

strains_CDE_rear = strains_rear(1:3, :);  % Group CDE (C: X-axis, D: 45°, E: Y-axis)
strains_FGH_rear = strains_rear(4:6, :);  % Group FGH (F: X-axis, G: 45°, H: Y-axis)
strains_IJK_rear = strains_rear(7:9, :);  % Group IJK (I: X-axis, J: 45°, K: Y-axis)

% Step 3: Calculate principal strains and shear strains for each rosette
% Function to calculate principal strains
function [ex, ey, gxy] = calculate_principal_strains(ea, eb, ec)
    theta = 45;  % Angle between rosette gauges
    ex = (ea + ec) / 2 + ((ea - ec) / 2) * cos(2 * deg2rad(theta)) + eb * sin(2 * deg2rad(theta));
    ey = (ea + ec) / 2 - ((ea - ec) / 2) * cos(2 * deg2rad(theta)) - eb * sin(2 * deg2rad(theta));
    gxy = - (ea - ec) * sin(2 * deg2rad(theta)) + eb * cos(2 * deg2rad(theta));
end

% http://www.engineeringcorecourses.com/solidmechanics1/C8-strain-transformation/C8.1-equations-of-strain-transformation/theory/

% Compute strains for each group
[ex_CDE_front, ey_CDE_front, gxy_CDE_front] = calculate_principal_strains(strains_CDE_front(1, :), strains_CDE_front(2, :), strains_CDE_front(3, :));
[ex_FGH_front, ey_FGH_front, gxy_FGH_front] = calculate_principal_strains(strains_FGH_front(1, :), strains_FGH_front(2, :), strains_FGH_front(3, :));
[ex_IJK_front, ey_IJK_front, gxy_IJK_front] = calculate_principal_strains(strains_IJK_front(1, :), strains_IJK_front(2, :), strains_IJK_front(3, :));

[ex_CDE_rear, ey_CDE_rear, gxy_CDE_rear] = calculate_principal_strains(strains_CDE_rear(1, :), strains_CDE_rear(2, :), strains_CDE_rear(3, :));
[ex_FGH_rear, ey_FGH_rear, gxy_FGH_rear] = calculate_principal_strains(strains_FGH_rear(1, :), strains_FGH_rear(2, :), strains_FGH_rear(3, :));
[ex_IJK_rear, ey_IJK_rear, gxy_IJK_rear] = calculate_principal_strains(strains_IJK_rear(1, :), strains_IJK_rear(2, :), strains_IJK_rear(3, :));

% Step 4: Convert principal strains to normal strains (epsilon_n)
theta = 45;  % Angle for normal strain calculation (in degrees)
theta_rad = deg2rad(theta);  % Convert angle to radians

% Function to calculate normal strain (epsilon_n)
function epsilon_n = calculate_normal_strain(e_x, e_y, g_xy, theta_rad)
    epsilon_n = e_x * cos(theta_rad)^2 + e_y * sin(theta_rad)^2 + g_xy * sin(theta_rad) * cos(theta_rad);
end

epsilon_en_CDE_front = calculate_normal_strain(ex_CDE_front, ey_CDE_front, gxy_CDE_front, theta_rad);
epsilon_en_FGH_front = calculate_normal_strain(ex_FGH_front, ey_FGH_front, gxy_FGH_front, theta_rad);
epsilon_en_IJK_front = calculate_normal_strain(ex_IJK_front, ey_IJK_front, gxy_IJK_front, theta_rad);

epsilon_en_CDE_rear = calculate_normal_strain(ex_CDE_rear, ey_CDE_rear, gxy_CDE_rear, theta_rad);
epsilon_en_FGH_rear = calculate_normal_strain(ex_FGH_rear, ey_FGH_rear, gxy_FGH_rear, theta_rad);
epsilon_en_IJK_rear = calculate_normal_strain(ex_IJK_rear, ey_IJK_rear, gxy_IJK_rear, theta_rad);


% Interpolate rear spar normal strain to match front spar load cases
epsilon_n_rear_interp_CDE = interp1(load_cases_rear, epsilon_en_CDE_rear, load_cases_front, 'linear');
epsilon_n_rear_interp_FGH = interp1(load_cases_rear, epsilon_en_FGH_rear, load_cases_front, 'linear');
epsilon_n_rear_interp_IJK = interp1(load_cases_rear, epsilon_en_IJK_rear, load_cases_front, 'linear');


% Step 5: Compute strain difference (front - rear)
strain_diff_CDE = epsilon_en_CDE_front - epsilon_n_rear_interp_CDE;
strain_diff_FGH = epsilon_en_FGH_front - epsilon_n_rear_interp_FGH;
strain_diff_IJK = epsilon_en_IJK_front - epsilon_n_rear_interp_IJK;

% Step 6: Compute torque (P * L) for each load case
torque = load_cases_front * L;  % Torque = Load * Moment Arm

% % Step 7: Plot results for both load cases and torque
% figure;
% subplot(2, 1, 1);
% hold on;
% plot(load_cases_front, epsilon_n_front, 'b-o', 'DisplayName', 'Front Spar (Normal Strain)');
% plot(load_cases_front, epsilon_n_rear_interp, 'r-o', 'DisplayName', 'Rear Spar (Normal Strain)');
% xlabel('Load (lbs)');
% ylabel('Normal Strain');
% title('Normal Strain vs Load for Front and Rear Spars');
% legend('Location', 'best');
% grid on;
% hold off;

% Plot Normal Strain Difference vs Torque (Load * Moment Arm)
subplot(3, 1, 1);
plot(torque, strain_diff_CDE, 'k-o', 'DisplayName', 'Normal Strain Difference');
xlabel('Torque (N * mm)');
ylabel('Normal Strain Difference');
title('Normal Strain Difference (Front - Rear) CDE vs Torque');
legend('Location', 'best');
grid on;

subplot(3, 1, 2);
plot(torque, strain_diff_FGH, 'k-o', 'DisplayName', 'Normal Strain Difference');
xlabel('Torque (N * mm)');
ylabel('Normal Strain Difference');
title('Normal Strain Difference (Front - Rear) FGH vs Torque');
legend('Location', 'best');
grid on;

subplot(3, 1, 3);
plot(torque, strain_diff_IJK, 'k-o', 'DisplayName', 'Normal Strain Difference');
xlabel('Torque (N * mm)');
ylabel('Normal Strain Difference');
title('Normal Strain Difference (Front - Rear) IJK vs Torque');
legend('Location', 'best');
grid on;

% Step 8: Plot Principal Strains for Each Rosette Group
% figure;
% hold on;
% plot(load_cases_front, ex_CDE, 'r-o', 'DisplayName', 'ex (CDE)');
% plot(load_cases_front, ey_CDE, 'b-o', 'DisplayName', 'ey (CDE)');
% plot(load_cases_front, gxy_CDE, 'g-o', 'DisplayName', 'gxy (CDE)');
% 
% plot(load_cases_front, ex_FGH, 'r--o', 'DisplayName', 'ex (FGH)');
% plot(load_cases_front, ey_FGH, 'b--o', 'DisplayName', 'ey (FGH)');
% plot(load_cases_front, gxy_FGH, 'g--o', 'DisplayName', 'gxy (FGH)');
% 
% plot(load_cases_front, ex_IJK, 'r-.o', 'DisplayName', 'ex (IJK)');
% plot(load_cases_front, ey_IJK, 'b-.o', 'DisplayName', 'ey (IJK)');
% plot(load_cases_front, gxy_IJK, 'g-.o', 'DisplayName', 'gxy (IJK)');
% 
% xlabel('Load (lbs)');
% ylabel('Strain');
% title('Principal and Shear Strains for Each Rosette Group');
% legend('Location', 'best');
% grid on;
% hold off;

% Save results for further analysis
% save('torque_normal_strain_data.mat', 'load_cases_front', 'strain_diff', 'torque');


% confused as the rosettes arent placed where we expect them to be