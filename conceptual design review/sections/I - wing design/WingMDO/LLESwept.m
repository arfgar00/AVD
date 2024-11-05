%% 
mywing = WingGeometry();
mywing.cr = 10;
mywing.ck = 8;
mywing.ct = 4;
mywing.s = 28.45/2;
mywing.Lambdain50 = 20*pi/180;
mywing.Lambdaout50 = 30*pi/180;
mywing.yk = 3;
mywing = mywing.calcSref();
mywing.N = 31;
figure(1)
clf;
mywing.plotWing()

airfoil = Airfoil();
airfoil = airfoil.readPolar("xf-sc21010-il-1000000.csv");
airfoil = airfoil.readShape("sc21010.dat.txt");
airfoil = airfoil.interpShape(9);
myairfoil = airfoil

c = mywing.cbar;
cruise = AirCondition()
cruise.M = 0.8;
cruise.h = convlength(39000, 'ft','m');
cruise = cruise.init(mywing.cbar)
myAirCondition = cruise
%% 
N = 50;
Damping = 2e-2

%% 
% Step 1: Define wing and flight parameters for a Boeing 777 
b = 2*mywing.s;
% Define total number of points


% Define the wing semi-span
s = mywing.s;  % Semi-span in meters

% Integration grid (y_i), from -s to s
y = linspace(-s, s, N)';  % Column vector
% Evaluation grid (y_n), offset from y_i to avoid overlap

% Adjust for the sweep angle by modifying the effective chord lengths
for i = 1:length(y)
    cn(i) = mywing.c_at_y(abs(y(i)));
    Lambda(i) = mywing.Lambdax_c(0.5,abs(y(i)));
end

% Display or use the chord distribution as needed
alpha = deg2rad(6);   % Angle of attack in radians
U_inf = myAirCondition.V;          % Freestream velocity in m/s (typical cruise speed)
rho = myAirCondition.rho;          % Air density in kg/m^3
converr = 1e-5;       % Convergence criteria

% Step 3: Initial guess for elliptical circulation distribution
Gamma_old = zeros(N,1); % Initialize Gamma_old to zeros
for i = 1:N
    if abs(y(i)) <= b/2
        Gamma_old(i) = sqrt(1 - (2*y(i)/b)^2); % Assign elliptical lift distribution
    else
        Gamma_old(i) = 0; % Set to zero outside the span
    end
end
Gamma_old = Gamma_old*800;
figure(2)
clf;
plot(y,Gamma_old)
hold on
plot(y, gradient(Gamma_old,y))
% Step 4 to 8: Iterative convergence loop
converged = false;
iteration = 0; % Counter to keep track of iteration
while converged == false
    iteration = iteration + 1;
    
    % Compute the derivative of Gamma with respect to y
    dGamma_dy = gradient(Gamma_old, y);  % Result is (Ny x 1)
    dGamma_dy = dGamma_dy(:);  % Converts dGamma_dy to a column vector (N x 1)
    % Integrate to find alpha induced(y)
    % Define a small exclusion distance epsilon
    epsilon = 2.22e-16; % Adjust this value based on the scale of y and the accuracy required
    for i = 1:length(y)
        % Define the integrand, excluding values where |y(i) - y| < epsilon
        integrand = dGamma_dy ./ (y(i) - y);
        
        % Apply the exclusion condition
        mask = abs(y - y(i)) > epsilon;
        integrand(~mask) = 0; % Set the integrand to zero where |y(i) - y| < epsilon
        
        % Perform the integration with trapz on the masked integrand
        alpha_induced(i) = 1 / (4 * pi * U_inf) * trapz(y(mask), integrand(mask));
    end
    
    %figure(2)
    %clf;
    %plot(y,alpha_induced)
    alphae = alpha - alpha_induced;
    for i = 1:length(alphae)
        CLp(i) = myairfoil.interpPolar(alphae(i));
    end
    Gamma = 1/2*U_inf.*cn'.*CLp'.*Lambda';
    %Gamma = 1/2*U_inf.*cn'.*CLp';

    figure(3)
    clf;
    plot(y,Gamma)
    hold on
    plot(y,Gamma_old)
    

    % Step 9: Check convergence
    diff = abs(Gamma_old - Gamma) ./ (abs(Gamma) ); 
    max_diff = max(diff); % Maximum relative difference for convergence check
    %fprintf('Iteration %d: Max relative difference = %.20f\n', iteration, max_diff);
    
    if max_diff < converr
        converged = true;
        disp('Process is complete! Converged circulation distribution found.');
        break;
    end

    if iteration > 100000
        break;
    end

    Gamma_old = Gamma_old + Damping * (Gamma - Gamma_old);
end
CL = 2/(U_inf * mywing.SREF) * trapz(y,Gamma);
CDi = 2/(U_inf * mywing.SREF) * trapz(y,Gamma.*alpha_induced);


figure(5)
clf;
plot(y,rho*U_inf*Gamma_old)
xlabel("y")
ylabel("Lift")