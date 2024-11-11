function [CL, CDi, CLy] = LLEKuchemann(mywing, myairfoil, myAirCondition, Damping)
    z = mywing.stripy; % Center of strips, the variable of integration
    
    % Initial guess 
    Gammaz = zeros(size(z)); % Initialize Gamma_old to zeros
    for i = 1:length(Gammaz)
        if abs(z(i)) <= mywing.s
            Gammaz(i) = 1*sqrt(1 - (2*z(i)/mywing.b)^2); % Assign elliptical lift distribution
        else
            Gammaz(i) = 0; % Set to zero outside the span
        end
    end
    %Gammaz_uncorrected = Gammaz;
    converged = false;
    epsilon = 2.22e-16;
    iteration = 0; % Counter to keep track of iteration
    converr = 1e-6;
    while converged == false
        iteration = iteration + 1;
        % Compute the derivative of Gamma with respect to y
        dGammaz_dz = gradient(Gammaz, z);  % 1xN

        for i = 1:length(z)
            % Define the integrand, excluding values where |y(i) - y| < epsilon
            integrand = dGammaz_dz ./ (z(i) -z);
            
            % Apply the exclusion condition
            mask = abs(z - z(i)) > epsilon;
            integrand(~mask) = 0; % Set the integrand to zero where |y(i) - y| < epsilon
            
            % Perform the integration with trapz on the masked integrand
            % Due to sweep, freestream velocity has become U_inf * cos(Lambda),
            % where Lambda is the sweep angle at quarter chords
            alpha_induced(i) = 1 / (4 * pi * myAirCondition.V * cos(mywing.Lambdax_c(0.25,abs(z(i))))) * trapz(z(mask), integrand(mask));
            %alpha_induced_uncorrected(i)= 1 / (4 * pi * myAirCondition.V) * trapz(z(mask), integrand(mask));
        end
        
        alphae = mywing.twistfun(z) - alpha_induced;
        %alphae_uncorrected = mywing.twistfun(z) - alpha_induced_uncorrected;
        for i = 1:length(z)
            % use modified lift curve slope for swept
            Lambdai = mywing.Lambdax_c(0.25,abs(z(i))); %Sweep angle for this index i
            M = myAirCondition.M; %free stream Mach number
            a0 = myairfoil.a0; %Original lift curve slope of airfoil.
            AR = mywing.AR;
            acomp = a0*cos(Lambdai) / (sqrt(1 - M^2 * cos(Lambdai)^2 + (a0*cos(Lambdai)/(pi*AR))^2) + a0*cos(Lambdai)/(pi*AR));
            CLp_Swept(i) = acomp * alphae(i) + myairfoil.b;
            %CLp_uncorrected(i) = myairfoil.interpPolar(alphae_uncorrected(i));
        end
        Gammaz_new = 1/2*myAirCondition.V.*mywing.cn.*CLp_Swept;
        %Gammaz_uncorrected_new =  1/2*myAirCondition.V.*mywing.cn.*CLp_uncorrected;
        
        Gammaz = Gammaz + Damping * (Gammaz_new - Gammaz);
        %Gammaz_uncorrected = Gammaz_uncorrected + Damping * (Gammaz_uncorrected_new - Gammaz_uncorrected);
        
        % figure(4)
        % clf;
        % plot(z,Gammaz)
        % hold on
        % plot(z,Gammaz_new)
        % hold on
        % plot(z,Gammaz_uncorrected,'--')
        % hold on
        % plot(z,Gammaz_uncorrected_new,'--')
        % legend("Gammaz","Gammaz new","Gammaz_uncorrected","Gammaz_uncorrected new")
        
    
        % Step 9: Check convergence
        diff = abs(Gammaz_new - Gammaz) ./ (abs(Gammaz) ); 
        max_diff = max(diff); % Maximum relative difference for convergence check
        %fprintf('Iteration %d: Max relative difference = %.20f\n', iteration, max_diff);
        
        if max_diff < converr
            converged = true;
            disp('Process is complete! Converged circulation distribution found.');
            break;
        end
    
        if iteration > 10000
            break;
        end
       
    end
    
    CL = 2/(myAirCondition.V * mywing.SREF) * trapz(z,Gammaz);
    CDi = 2/(myAirCondition.V * mywing.SREF) * trapz(z,Gammaz.*alpha_induced);
    % 
    % CL_uncorrected = 2/(myAirCondition.V * mywing.SREF) * trapz(z,Gammaz_uncorrected);
    % CDi_uncorrected = 2/(myAirCondition.V * mywing.SREF) * trapz(z,Gammaz_uncorrected.*alpha_induced_uncorrected);
    
    
    %L = AirCondition.rho*myAirCondition.V*Gammaz;
    Cldistribution = Gammaz./(1/2*mywing.SREF*myAirCondition.V);
    
    % figure(3)
    % clf;
    % plot(z, Cldistribution)
    CLy = [z; Cldistribution];
end