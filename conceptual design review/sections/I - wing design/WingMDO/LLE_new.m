function [CL,CDi,Cly] = LLE_new(mywing,myairfoil,myAirCondition,N,Damping)
    %linear twist

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
        Lambda(i) = mywing.Lambdax_c(0.25,abs(y(i))); %sweep at quarter quard
        alpha(i) = mywing.twist_max / mywing.s * abs(y(i));
    end
    
    % Display or use the chord distribution as needed
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
        Gamma = 1/2*U_inf.*cn'.*CLp'.*cos(Lambda');
    
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
    
        Gamma_old = Gamma_old + Damping * (Gamma - Gamma_old);
    end
    CL = 2/(U_inf * mywing.SREF) * trapz(y,Gamma);
    CDi = 2/(U_inf * mywing.SREF) * trapz(y,Gamma.*alpha_induced');
    
    L = rho*U_inf*Gamma;
    Cldistribution = L./(1/2*rho*mywing.SREF*U_inf^2);
    Cly = [y Cldistribution];
    Cly = Cly';
end