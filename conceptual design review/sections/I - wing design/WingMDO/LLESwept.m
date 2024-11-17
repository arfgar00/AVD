function [CL, CDi, CLy, CDprofile] = LLESwept(mywing, myairfoil, myAirCondition)
    global Damping tolerance bodyDiameter
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
    %No lift at fuselage
    %indices = (y >= -bodyDiameter/2) & (y <= bodyDiameter/2);
    %Gamma(indices) = 0;

    alpha_i = zeros(1, k);   % Initial guess for induced angle of attack
    %alpha_i_uncorrected = alpha_i;   % Initial guess for induced angle of attack
    
    
    % Iterative solver for the circulation distribution
    %tolerance = 1e-4; % Convergence tolerance
    error = inf;       % Pre-define error
    iteration = 0;     % Number of iteration
    while error > tolerance %|| error_uncorrected > tolerance
        Gamma_old = Gamma; % Define guessed circulation
        %Gamma_old_uncorrected = Gamma_uncorrected; 
        % Compute induced angle of attack (alpha_i) at each station
    
        for n = 1 : k
            Lambdai = mywing.Lambdax_c(0.25, y(n));
            integral = 0; % Pre-define integral
            %integral_uncorrected = 0;
            for j = 1 : k
                if j ~= n
                    % Compute dGamma/dy using central differencing scheme
                    if j == 1
                        dGammady = (Gamma(j + 1) - Gamma(j)) / dy;
                        %dGammady_uncorrected = (Gamma_uncorrected(j + 1) - Gamma_uncorrected(j)) / dy;
                    elseif j == k
                        dGammady = (Gamma(j) - Gamma(j - 1)) / dy;
                        %dGammady_uncorrected = (Gamma_uncorrected(j) - Gamma_uncorrected(j - 1)) / dy;
                    else
                        dGammady = (Gamma(j + 1) - Gamma(j - 1)) / (2 * dy);
                        %dGammady_uncorrected = (Gamma_uncorrected(j + 1) - Gamma_uncorrected(j - 1)) / (2 * dy);
                    end
                    
                    % Handle singularity by averaging adjacent sections
    
                    if abs(y(n) - y(j)) < 1e-6
                        if j == 1
                            integral = integral + dGammady / ((y(n) - y(j + 1)) / 2);
                            %integral_uncorrected = integral_uncorrected + dGammady_uncorrected / ((y(n) - y(j + 1)) / 2);
                        elseif j == k
                            integral = integral + dGammady / ((y(n) - y(j - 1)) / 2);
                            %integral_uncorrected = integral_uncorrected + dGammady_uncorrected / ((y(n) - y(j - 1)) / 2);
                        else
                            integral = integral + dGammady / ((y(n) - y(j - 1)) / 2 + (y(n) - y(j + 1)) / 2);
                            %integral_uncorrected = integra_uncorrectedl + dGammady_uncorrected / ((y(n) - y(j - 1)) / 2 + (y(n) - y(j + 1)) / 2);
    
                        end
    
                    else
                        integral = integral + dGammady / (y(n) - y(j));
                        %integral_uncorrected = integral_uncorrected + dGammady_uncorrected / (y(n) - y(j));
                    end
    
                end
                
            end
    
            %alpha_i_uncorrected(n) = (1 / (4 * pi * U_inf)) * dy * integral_uncorrected;
            alpha_i(n) = (1 / (4 * pi * U_inf * cos(Lambdai))) * dy * integral;
        end
        
        % Calculate effective angle of attack
        alpha_e = mywing.twistfun(y) - alpha_i;
        %alpha_e_uncorrected = mywing.twistfun(y) - alpha_i_uncorrected;
        % Update circulation distribution using sectional lift coefficient
        for n = 1:k
            Lambdai = mywing.Lambdax_c(0.25, y(n));
            c_n = mywing.cn(n);
            acomp = a0*cos(Lambdai) / (sqrt(1 - M^2 * cos(Lambdai)^2 + (a0*cos(Lambdai)/(pi*AR))^2) + a0*cos(Lambdai)/(pi*AR));
            acomp = real(acomp);
            %CL_n_uncorrected = a0 * alpha_e_uncorrected(n) + b;
            %STall is considered
            if (alpha_e(n) > myairfoil.alphaSPos) || (alpha_e(n) < myairfoil.alphaSNeg)
                CL_n(n) = 0;
            else
                CL_n(n) = acomp * alpha_e(n) + b;
            end
                
            %Gamma_uncorrected(n) = 0.5 * U_inf * c_n * CL_n_uncorrected;
            Gamma(n) =  0.5 * U_inf * cos(Lambdai) * c_n * CL_n(n);
        end
        
        % Apply boundary conditions again
        Gamma(1) = 0; % Circulation at the wing tip (y = -s) is zero
        Gamma(k) = 0; % Circulation at the wing tip (y = s) is zero
        %Gamma_uncorrected(1) = 0; % Circulation at the wing tip (y = -s) is zero
        %Gamma_uncorrected(k) = 0; % Circulation at the wing tip (y = s) is zero
        
        % Damping iteration
        Gamma = Gamma_old + Damping * (Gamma - Gamma_old);
        %Gamma_uncorrected = Gamma_old_uncorrected + Damping * (Gamma_uncorrected - Gamma_old_uncorrected);
       
        % Calculate error for convergence
        %error_uncorrected = max(abs(Gamma_uncorrected - Gamma_old_uncorrected));
        error = max(abs(Gamma - Gamma_old));
        iteration = iteration + 1; 
        %disp(["error = ", error])
        figure(3)
        clf;
        plot(y,Gamma)
        hold on
        plot(y,Gamma_old)
        xlabel("y")
        ylabel("Gamma")
        legend("Gamma","Gamma_old")
    end
    
    indices = (y >= -bodyDiameter/2) & (y <= bodyDiameter/2);
    CL_n(indices) = 0;
    
    
    % Calculate results using Simpson's rule
    CL = (2 / (U_inf * S_ref)) * (dy / 3) * (Gamma(1) + ...
        4 * sum(Gamma(2 : 2 : end - 1)) + 2 * sum(Gamma(3 : 2 : end - 2)) + Gamma(end));
    %CL_uncorrected = (2 / (U_inf * S_ref)) * (dy / 3) * (Gamma_uncorrected(1) + ...
        %4 * sum(Gamma_uncorrected(2 : 2 : end - 1)) + 2 * sum(Gamma_uncorrected(3 : 2 : end - 2)) + Gamma_uncorrected(end));
    CDi = (2 / (U_inf * S_ref)) * (dy / 3) * (Gamma(1) * alpha_i(1) + ...
        4 * sum(Gamma(2 : 2 : end - 1) .* alpha_i(2 : 2 : end - 1)) + ...
        2 * sum(Gamma(3 : 2 : end - 2) .* alpha_i(3 : 2 : end - 2)) + Gamma(end) * alpha_i(end));
    CDprofile = 0;
    for i = 1:length(y)
        CDprofile = CDprofile + interp1(myairfoil.alpha,myairfoil.Cd,alpha_e(i)).*mywing.Sc(i)./mywing.SREF;
    end
    %Cldistribution_uncorrected = Gamma_uncorrected./(1/2*mywing.SREF*myAirCondition.V);
    %disp([size(CL_n)])
    for i = 1:length(y)
        localChord(i) = mywing.c_at_y(y(i));
    end
    CL_n(1) = CL_n(2)*localChord(1)/localChord(2);
    CL_n(k) = CL_n(1);
    CLy = [y; CL_n.*localChord./mywing.cbar]; %2xN
    figure(5)
    clf;
    plot(y,CL_n.*localChord./mywing.cbar)
    xlabel("y")
    ylabel("Cl")
end