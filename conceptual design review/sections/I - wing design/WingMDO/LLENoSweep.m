%% 
clear all
mywing = WingGeometry();
mywing.cr = 1;
mywing.ck = 2/3;
mywing.ct = 1/3;
mywing.s = 4;
mywing.Lambdain50 = 0;
mywing.Lambdaout50 = 0;
mywing.yk = 2;
mywing.N = 501;
mywing.twist_max = deg2rad(5);
mywing = mywing.calcSref();
mywing = mywing.createStrips();
figure(1)
clf;
mywing.plotWing()
axis equal

figure(2)
clf;
myairfoil = Airfoil();
myairfoil = myairfoil.readPolar("xf-sc21010-il-1000000.csv");
myairfoil = myairfoil.readShape("sc21010.dat.txt");
myairfoil = myairfoil.interpShape(9);
myairfoil.plotPolar()


myAirCondition = AirCondition();
myAirCondition.M = 0.83;
myAirCondition.h = convlength(39000, 'ft','m');
myAirCondition = myAirCondition.init(mywing.cbar);
myAirCondition.V = 50;
myAirCondition = myAirCondition.calcM();
myAirCondition.rho = 1.225;

DFuselage = 6;
%% Method in ref26 and numerical LLE notes
D = 0.01;          % Damping factor
s = mywing.s;
alpha = @(y) deg2rad(5) * (1 - (y.^2 / s^2)); % Twist angle distribution

k = mywing.N;
U_inf = myAirCondition.V;
AR = mywing.AR;
S_ref = mywing.SREF;
c_avg = mywing.cbar;
y = mywing.stripy; % Center of strips, the variable of integration
dy = (2 * mywing.s) / mywing.N;       % Step size for y
c_values = mywing.cn;
%c_values = c(y); % Evaluate chord length at all spanwise locations

% Initial guess
Gamma = 0.1 * ones(1, k); % Initial guess for circulation distribution
alpha_i = zeros(1, k);   % Initial guess for induced angle of attack

% Set boundary condition
Gamma(1) = 0; % circulation at the wing tip (y = -s) is zero
Gamma(k) = 0; % circulation at the wing tip (y = s) is zero

% Iterative solver for the circulation distribution

tolerance = 1e-13; % Convergence tolerance
error = inf;       % Pre-define error
iteration = 0;     % Number of iteration
a0 = 6;
while error > tolerance

    Gamma_old = Gamma; % Define guessed circulation
    
    % Compute induced angle of attack (alpha_i) at each station

    for n = 1 : k

        integral = 0; % Pre-define integral

        for j = 1 : k

            if j ~= n

                % Compute dGamma/dy using central differencing scheme

                if j == 1

                    dGammady = (Gamma(j + 1) - Gamma(j)) / dy;

                elseif j == k

                    dGammady = (Gamma(j) - Gamma(j - 1)) / dy;

                else

                    dGammady = (Gamma(j + 1) - Gamma(j - 1)) / (2 * dy);

                end
                
                % Handle singularity by averaging adjacent sections

                if abs(y(n) - y(j)) < 1e-6

                    if j == 1

                        integral = integral + dGammady / ((y(n) - y(j + 1)) / 2);

                    elseif j == k

                        integral = integral + dGammady / ((y(n) - y(j - 1)) / 2);

                    else

                        integral = integral + dGammady / ((y(n) - y(j - 1)) / 2 + (y(n) - y(j + 1)) / 2);

                    end

                else

                    integral = integral + dGammady / (y(n) - y(j));

                end

            end
            
        end

        alpha_i(n) = (1 / (4 * pi * U_inf)) * dy * integral;
        
    end
    
    % Calculate effective angle of attack
    alpha_e = alpha(y) - alpha_i;
    
    % Update circulation distribution using sectional lift coefficient
    for n = 1:k

        c_n = mywing.cn(n);

        CL_n = a0 * alpha_e(n);

        Gamma(n) = 0.5 * U_inf * c_n * CL_n;

    end
    
    % Apply boundary conditions again
    Gamma(1) = 0; % Circulation at the wing tip (y = -s) is zero
    Gamma(k) = 0; % Circulation at the wing tip (y = s) is zero
    
    % Damping iteration
    Gamma = Gamma_old + D * (Gamma - Gamma_old);
    
    % Calculate error for convergence
    error = max(abs(Gamma - Gamma_old));
    iteration = iteration + 1;
    
end

% Calculate results using Simpson's rule
CL = (2 / (U_inf * S_ref)) * (dy / 3) * (Gamma(1) + ...
    4 * sum(Gamma(2 : 2 : end - 1)) + 2 * sum(Gamma(3 : 2 : end - 2)) + Gamma(end));
CDi = (2 / (U_inf * S_ref)) * (dy / 3) * (Gamma(1) * alpha_i(1) + ...
    4 * sum(Gamma(2 : 2 : end - 1) .* alpha_i(2 : 2 : end - 1)) + ...
    2 * sum(Gamma(3 : 2 : end - 2) .* alpha_i(3 : 2 : end - 2)) + Gamma(end) * alpha_i(end));

% Display results
fprintf('Converged after %d iterations\n', iteration);

figure
% Plot wing shape, lift distribution, and induced angle distribution

% Plot wing shape
subplot(3, 1, 1);
plot(y, c_values, 'b', 'LineWidth', 2);
xlabel('Spanwise Location (y)');
ylabel('Chord Length (c)');
title('Wing Shape');
grid on;

% Plot lift distribution along the span
subplot(3, 1, 2);
plot(y, Gamma, 'r', 'LineWidth', 2);
xlabel('Spanwise Location (y)');
ylabel('Circulation (Gamma)');
title('Lift Distribution Along the Span');
grid on;

% Plot induced angle distribution along the span
subplot(3, 1, 3);
plot(y, rad2deg(alpha_i), 'g', 'LineWidth', 2);
xlabel('Spanwise Location (y)');
ylabel('Induced Angle of Attack (degrees)');
title('Induced Angle Distribution Along the Span');
grid on;