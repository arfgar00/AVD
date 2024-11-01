function [C_L,CDi] = LLE(mywing,myairfoil,myAirCondition,N,D)
        % Step 1: Define wing and flight parameters for a Boeing 777 
        b = 2*mywing.s;
        y = linspace(-mywing.s, mywing.s, N);  % Spanwise positions from -b/2 to b/2
        % Adjust for the sweep angle by modifying the effective chord lengths
        for i = 1:length(y)
            chord(i) = mywing.c_at_y(abs(y(i)));
        end
        
        % Display or use the chord distribution as needed
        alpha = deg2rad(12.81);   % Angle of attack in radians
        V_inf = myAirCondition.V;          % Freestream velocity in m/s (typical cruise speed)
        rho = myAirCondition.rho;          % Air density in kg/m^3
        epsilon = 1e-13;       % Convergence criteria
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
            %% Important change happened here
            [Cl,Cd] = myairfoil.interpPolar(alpha_eff);
            %a0 = 6.105;
            %Cl = a0 * alpha_eff;
        
            % Step 8: Calculate new circulation using Kutta-Joukowski theorem
            Gamma_new = Cl .* chord .* V_inf / 2;
        
            % Step 9: Check convergence
            diff = abs(Gamma_old - Gamma_new) ./ (abs(Gamma_new) ); 
            max_diff = max(diff); % Maximum relative difference for convergence check
            %fprintf('Iteration %d: Max relative difference = %.20f\n', iteration, max_diff);
            
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
        fprintf('Iteration %d: Max relative difference = %.20f\n', iteration, max_diff);
        % Step 11: Calculate total circulation using Simpson's rule
        total_circulation = (dy / 3) * (Gamma_old(1) + 4 * sum(Gamma_old(2:2:end)) + 2 * sum(Gamma_old(3:2:end))+ Gamma_old(end));
        total_circulation_alpha = (dy / 3) * (Gamma_old(1).*alpha_induced(1) + 4 * sum(Gamma_old(2:2:end).*alpha_induced(2:2:end)) + 2 * sum(Gamma_old(3:2:end).*alpha_induced(3:2:end))+ Gamma_old(end).*alpha_induced(end));
        

        % Step 12: Calculate reference area (S)
        S = (b) * mean(chord); % Reference area based on average chord length
        
        % Step 13: Calculate sectional lift coefficient C_L
        C_L = (2 / (V_inf * S)) * total_circulation;
        
        % Step 14: Plot the final circulation distribution across the span
        figure;
        clf;
        plot(y, Gamma_old, 'LineWidth', 1.5); % Plot circulation
        xlabel('Spanwise Position (y)');
        ylabel('Circulation (\Gamma)');
        grid on;
        
        CDi = (2 / (V_inf * S)) * total_circulation_alpha;
        
        % Step 17: Output the total lift coefficient C_L
        fprintf('Total Lift Coefficient (C_L): %.6f\n', C_L);
        fprintf('Total Induced Drag (CDi): %.6f\n', CDi);
end