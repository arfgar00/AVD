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
mywing.twistfun = @(y) deg2rad(5) * (1 - (y.^2 / mywing.s^2));
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
myairfoil.a0 = 6;
myairfoil.b = 0;


myAirCondition = AirCondition();
myAirCondition.M = 0.83;
myAirCondition.h = convlength(39000, 'ft','m');
myAirCondition = myAirCondition.init(mywing.cbar);
myAirCondition.V = 50;
myAirCondition = myAirCondition.calcM();
myAirCondition.rho = 1.225;

DFuselage = 6;
D = 0.01;          % Damping factor

Diameter = 2;
%% Method in ref26 and numerical LLE notes

s = mywing.s;
k = mywing.N;
S_ref = mywing.SREF;
c_avg = mywing.cbar;
y = mywing.stripy; % Center of strips, the variable of integration
dy = (2 * mywing.s) / mywing.N;       % Step size for y
c_values = mywing.cn;
AR = mywing.AR;

U_inf = myAirCondition.V;
M = myAirCondition.M;

a0 = myairfoil.a0;
b = myairfoil.b;

% Initial guess
Gamma = 0.1 * ones(1, k); % Initial guess for circulation distribution
% Set boundary condition
Gamma(1) = 0; % circulation at the wing tip (y = -s) is zero
Gamma(k) = 0; % circulation at the wing tip (y = s) is zero
Gamma_uncorrected = Gamma; % Initial guess for circulation distribution

alpha_i = zeros(1, k);   % Initial guess for induced angle of attack
alpha_i_uncorrected = alpha_i;   % Initial guess for induced angle of attack


% Iterative solver for the circulation distribution
tolerance = 1e-13; % Convergence tolerance
error = inf;       % Pre-define error
iteration = 0;     % Number of iteration
while error > tolerance || error_uncorrected > tolerance
    Gamma_old = Gamma; % Define guessed circulation
    Gamma_old_uncorrected = Gamma_uncorrected; 
    % Compute induced angle of attack (alpha_i) at each station

    for n = 1 : k
        Lambdai = mywing.Lambdax_c(0.25, y(n));
        integral = 0; % Pre-define integral
        integral_uncorrected = 0;
        for j = 1 : k
            if j ~= n
                % Compute dGamma/dy using central differencing scheme
                if j == 1
                    dGammady = (Gamma(j + 1) - Gamma(j)) / dy;
                    dGammady_uncorrected = (Gamma_uncorrected(j + 1) - Gamma_uncorrected(j)) / dy;
                elseif j == k
                    dGammady = (Gamma(j) - Gamma(j - 1)) / dy;
                    dGammady_uncorrected = (Gamma_uncorrected(j) - Gamma_uncorrected(j - 1)) / dy;
                else
                    dGammady = (Gamma(j + 1) - Gamma(j - 1)) / (2 * dy);
                    dGammady_uncorrected = (Gamma_uncorrected(j + 1) - Gamma_uncorrected(j - 1)) / (2 * dy);
                end
                
                % Handle singularity by averaging adjacent sections

                if abs(y(n) - y(j)) < 1e-6
                    if j == 1
                        integral = integral + dGammady / ((y(n) - y(j + 1)) / 2);
                        integral_uncorrected = integral_uncorrected + dGammady_uncorrected / ((y(n) - y(j + 1)) / 2);
                    elseif j == k
                        integral = integral + dGammady / ((y(n) - y(j - 1)) / 2);
                        integral_uncorrected = integral_uncorrected + dGammady_uncorrected / ((y(n) - y(j - 1)) / 2);
                    else
                        integral = integral + dGammady / ((y(n) - y(j - 1)) / 2 + (y(n) - y(j + 1)) / 2);
                        integral_uncorrected = integra_uncorrectedl + dGammady_uncorrected / ((y(n) - y(j - 1)) / 2 + (y(n) - y(j + 1)) / 2);

                    end

                else
                    integral = integral + dGammady / (y(n) - y(j));
                    integral_uncorrected = integral_uncorrected + dGammady_uncorrected / (y(n) - y(j));
                end

            end
            
        end

        alpha_i_uncorrected(n) = (1 / (4 * pi * U_inf)) * dy * integral_uncorrected;
        alpha_i(n) = (1 / (4 * pi * U_inf * cos(Lambdai))) * dy * integral;
    end
    
    % Calculate effective angle of attack
    alpha_e = mywing.twistfun(y) - alpha_i;
    alpha_e_uncorrected = mywing.twistfun(y) - alpha_i_uncorrected;
    % Update circulation distribution using sectional lift coefficient
    for n = 1:k
        Lambdai = mywing.Lambdax_c(0.25, y(n));
        c_n = mywing.cn(n);
        acomp = a0*cos(Lambdai) / (sqrt(1 - M^2 * cos(Lambdai)^2 + (a0*cos(Lambdai)/(pi*AR))^2) + a0*cos(Lambdai)/(pi*AR));
        CL_n_uncorrected = a0 * alpha_e_uncorrected(n) + b;
        CL_n = acomp * alpha_e(n) + b;

        Gamma_uncorrected(n) = 0.5 * U_inf * c_n * CL_n_uncorrected;
        Gamma(n) =  0.5 * U_inf * cos(Lambdai) * c_n * CL_n;
    end
    
    % Apply boundary conditions again
    Gamma(1) = 0; % Circulation at the wing tip (y = -s) is zero
    Gamma(k) = 0; % Circulation at the wing tip (y = s) is zero
    %indices = (y >= -Diameter/2) & (y <= Diameter/2);
    %Gamma(indices) = 0;
    Gamma_uncorrected(1) = 0; % Circulation at the wing tip (y = -s) is zero
    Gamma_uncorrected(k) = 0; % Circulation at the wing tip (y = s) is zero
    
    % Damping iteration
    Gamma = Gamma_old + D * (Gamma - Gamma_old);
    Gamma_uncorrected = Gamma_old_uncorrected + D * (Gamma_uncorrected - Gamma_old_uncorrected);
   
    % Calculate error for convergence
    error_uncorrected = max(abs(Gamma_uncorrected - Gamma_old_uncorrected));
    error = max(abs(Gamma - Gamma_old));
    iteration = iteration + 1; 
end

% Calculate results using Simpson's rule
CL = (2 / (U_inf * S_ref)) * (dy / 3) * (Gamma(1) + ...
    4 * sum(Gamma(2 : 2 : end - 1)) + 2 * sum(Gamma(3 : 2 : end - 2)) + Gamma(end));
CL_uncorrected = (2 / (U_inf * S_ref)) * (dy / 3) * (Gamma_uncorrected(1) + ...
    4 * sum(Gamma_uncorrected(2 : 2 : end - 1)) + 2 * sum(Gamma_uncorrected(3 : 2 : end - 2)) + Gamma_uncorrected(end));
CDi = (2 / (U_inf * S_ref)) * (dy / 3) * (Gamma(1) * alpha_i(1) + ...
    4 * sum(Gamma(2 : 2 : end - 1) .* alpha_i(2 : 2 : end - 1)) + ...
    2 * sum(Gamma(3 : 2 : end - 2) .* alpha_i(3 : 2 : end - 2)) + Gamma(end) * alpha_i(end));


Cldistribution_uncorrected = Gamma_uncorrected./(1/2*mywing.SREF*myAirCondition.V);

Cldistribution = Gamma./(1/2*mywing.SREF*myAirCondition.V);

% Display results
fprintf('Converged after %d iterations\n', iteration);

figure(3)
clf;
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

figure(4)
clf;
plot(y,Cldistribution)
hold on
plot(y, Cldistribution_uncorrected,'--')