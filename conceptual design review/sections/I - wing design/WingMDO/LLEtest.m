%% 
wing = WingGeometry();
wing.cr = 5;
wing.ck = 5;
wing.ct = 5;
wing.s = 28.45/2;
wing.Lambdain50 = 0*pi/180;
wing.Lambdaout50 = 0*pi/180;
wing.yk = 5.69;
wing = wing.calcSref();
wing.N = 31;
mywing = wing;

airfoil = Airfoil();
airfoil = airfoil.readPolar("xf-sc21010-il-1000000.csv");
airfoil = airfoil.readShape("sc21010.dat.txt");
airfoil = airfoil.interpShape(9);
myairfoil = airfoil

c = wing.cbar;
cruise = AirCondition()
cruise.M = 0.8;
cruise.h = convlength(39000, 'ft','m');
cruise = cruise.init(wing.cbar)
myAirCondition = cruise
%% 
N = 20;
Damping = 8e-2

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
end

% Display or use the chord distribution as needed
alpha = deg2rad(8);   % Angle of attack in radians
U_inf = myAirCondition.V;          % Freestream velocity in m/s (typical cruise speed)
rho = myAirCondition.rho;          % Air density in kg/m^3
epsilon = 1e-13;       % Convergence criteria

% Step 3: Initial guess for elliptical circulation distribution
Gamma_old = zeros(1, N); % Initialize Gamma_old to zeros
for i = 1:N
    if abs(y(i)) <= b/2
        Gamma_old(i) = sqrt(1 - (2*y(i)/b)^2); % Assign elliptical lift distribution
    else
        Gamma_old(i) = 0; % Set to zero outside the span
    end
end
Gamma_old = Gamma_old.*1000
figure(1)
clf;
plot(y,Gamma_old)
hold on
plot(y, gradient(Gamma_old,y))
% Step 4 to 8: Iterative convergence loop
converged = false;
iteration = 0; % Counter to keep track of iteration
while iteration <= 10000
    iteration = iteration + 1;
    
    % Compute the derivative of Gamma with respect to y
    dGamma_dy = gradient(Gamma_old, y);  % Result is (Ny x 1)
    dGamma_dy = dGamma_dy(:);  % Converts dGamma_dy to a column vector (N x 1)
    % Step 2: Create the difference matrix D
    % D(i, j) = y_i - y_j
    [Y_i, Y_j] = meshgrid(y, y);   % Both are (N x N)
    D = Y_i - Y_j;                 % (N x N)
    
    % Step 3: Handle the singularity at D = 0
    % Set diagonal elements to NaN or a small value to avoid division by zero
    D(logical(eye(N))) = NaN;  % Alternatively, D(logical(eye(N))) = eps;
    
    % Step 4: Replicate dGamma_dy for matrix operations
    % We need dGamma_dy(j) for each (i, j)
    dGamma_dy_mat = ones(N, 1) * dGamma_dy';
    
    % Step 5: Compute the integrand
    integrand = dGamma_dy_mat ./ D;  % Element-wise division, (N x N)
    
    % Step 6: Handle the singularity in the integrand
    integrand(isnan(integrand)) = 0;  % Set NaN values to zero
    
    % Step 7: Perform numerical integration over y_j for each y_i
    % Integrate along the columns (dimension 2)
    alpha_induced = (1 / (4 * pi * U_inf)) * trapz(y, integrand, 2);  % (N x 1)
    %figure(2)
    %clf;
    %plot(y,alpha_induced)
    alphae = alpha - alpha_induced;
    CLp = myairfoil.interpPolar(alphae);
    Gamma = 1/2*U_inf.*cn'.*CLp;
    Gamma_old = Gamma_old(:);
    figure(3)
    clf;
    plot(y,Gamma)
    hold on
    plot(y,Gamma_old)
    Gamma_old = Gamma_old + Damping * (Gamma - Gamma_old);
end