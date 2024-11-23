function [a,b] = wingLiftCurveSlope(wing,airfoil,h,M,a1,a2,N,plotgraph,plotstall)
    %inputs:
    % wing: a WingGeometry object
    % airfoil: a Airfoil object
    % h: altitude in meter
    % M: Mach number
    % a1: Start of angle of attack in radians
    % a2: End of angle of attack in radians
    % N: Number of interpolations, eg. 10
    % plotgraph: if set to true will plot CLvs.alpha
    % plotstall: if set to true will plot stall visualization
    %%
    function [CLlist, CDilist, CDprofilelist, stallylist] = LLESwept_aseq(mywing, myairfoil, myAirCondition,Damping, tolerance, bodyDiameter, aseq)
        % compute CL at angle of attack specified in aseq, list of alpha
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
        stallAlpha = 0;
        stallylist = zeros(mywing.N,length(aseq));
        for Idx = 1:length(aseq)
            alpha = aseq(Idx);
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
                alpha_e = alpha + mywing.twistfun(y) - alpha_i;
    
                % Update circulation distribution using sectional lift coefficient
                for n = 1:k
                    Lambdai = mywing.Lambdax_c(0.25, y(n));
                    c_n = mywing.cn(n);
                    acomp = a0*cos(Lambdai) / (sqrt(1 - M^2 * cos(Lambdai)^2 + (a0*cos(Lambdai)/(pi*AR))^2) + a0*cos(Lambdai)/(pi*AR));
                    acomp = real(acomp);
                    %STall is considered
                    if (alpha_e(n) > myairfoil.alphaSPos) || (alpha_e(n) < myairfoil.alphaSNeg)
                        CL_n(n) = 0;
                        if abs(mywing.stripy(n)) > bodyDiameter/2 %% stall of an actual section not inside fuselage
                            if stallAlpha == 0 && (n~=1) && (n~=mywing.N)
                                stallAlpha = alpha;
                                disp(["Stall at alpha = ", alpha.*180/pi])
                                disp(["Stall position y = ", mywing.stripy(n)])
                                disp(["alpha_e", alpha_e(n)*180/pi])
                            end
                            stallylist(n,Idx) = 1;
                            
                        end
                    else
                        CL_n(n) = acomp * alpha_e(n) + b;
                    end
                        
                    %Gamma_uncorrected(n) = 0.5 * U_inf * c_n * CL_n_uncorrected;
                    Gamma(n) =  0.5 * U_inf * cos(Lambdai) * c_n * CL_n(n);
                end
                
                % Apply boundary conditions again
                Gamma(1) = 0; % Circulation at the wing tip (y = -s) is zero
                Gamma(k) = 0; % Circulation at the wing tip (y = s) is zero
                % Damping iteration
                Gamma = Gamma_old + Damping * (Gamma - Gamma_old);
               
                % Calculate error for convergence
                error = max(abs(Gamma - Gamma_old));
                iteration = iteration + 1; 
            end
       
        %% End of LLE function
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
            if abs(y) > bodyDiameter/2
                CDprofile = CDprofile + interp1(myairfoil.alpha,myairfoil.Cd,alpha_e(i)).*mywing.Sc(i)./mywing.SREF;
            end
        end

        CLlist(Idx) = CL;
        CDilist(Idx) = CDi;
        CDprofilelist(Idx) = CDprofile;
    end
end
    aseq = linspace(a1,a2,N);
    Damping = 0.04;
    tolerance = 1e-2;
    bodyDiameter = 6.38;
    cruise1 = AirCondition();
    cruise1.M = M;
    cruise1.h = h;
    cruise1 = cruise1.init(wing.cbar);
    [CLlist_1, CDilist_1, CDprofilelist_1, stallyList] = LLESwept_aseq(wing, airfoil, cruise1 ,Damping, tolerance, bodyDiameter, aseq)
    
    [max_CL, max_idx] = max(CLlist_1);
    
    % Find the corresponding aseq
    aseq_at_max_CL = aseq(max_idx);
    % Interpolation range
    CL_range = [CLlist_1(1), max(CLlist_1)];
    aseq_range = [aseq(1), aseq_at_max_CL]; % Corresponding 'aseq' values
    
    % Fit a straight line using polyfit
    p = polyfit([aseq_range(1), aseq_range(2)], [CL_range(1), CL_range(2)], 1);
    
    % Generate interpolated values within the range
    aseq_interp = linspace(aseq_range(1), aseq_range(2), 100); % Interpolated x-values
    CL_interp = polyval(p, aseq_interp); % Evaluate the straight-line fit
    a = p(1);
    b = p(2);
    if plotgraph == true
        figure
        clf;
        scatter(aseq*180/pi, CLlist_1)
        % Plot the interpolated line
        hold on
        plot(aseq_interp*180/pi, CL_interp, '-r', 'LineWidth', 1.5);
        ylabel("CL")
        xlabel("alpha [degrees]")
        %title(getVariableName(cruise1))
    end
    if plotstall == true
        figure clf;
    
        for i =1:length(aseq)
            figure
            clf;
            wing_opt.plotWing()
            hold on
            wing_opt.plotStallOnStrips(stallyList(:,i))
            ylim([0 65/2])
            title(["alpha = ", aseq(i)*180/pi])
        end
    end
end