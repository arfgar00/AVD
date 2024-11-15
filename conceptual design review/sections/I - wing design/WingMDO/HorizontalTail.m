classdef HorizontalTail
    properties (Dependent)
        TailDimensions                     % Structure containing B_ht, C_ht_tip, C_ht_root, and MAC
        S_ht                               % Planform area (S_ht) in square meters
        L_ht                               % Tail arm (L_ht) in meters
        incidence_i_ht                     % Incidence angle (i_ht) in rad
        CM_owf                             % Wing-Fuselage pitching moment coefficient (CM_owf) in unitless
        CL_h                               % Coefficient of Lift of the horizontal tailplane (CL_h) in unitless
        CL_c                               % Coefficient of Lift at cruise (CL_c) in unitless
        CL_alphaH                          % Horizontal lift curve slope (C_Lα) in unitless 
        epsilon                            % Downwash angle (ε) in rad
        alphaH_eff                         % Effective angle of attack of the tailplane (αh) in degrees
        CLH_c                              % Horizontal tail lift coefficient at cruise (CLH_c)
    end

    properties (SetAccess = public)
        % Aircraft Geometry & Aerodynamics
        C_bar = 1;                         % Wing Mean Aerodynamic Chord (C_bar) in m
        S_w = 1;                           % Wing Reference Area (S) in m^2
        Df = 1;                            % Fuselage Diameter (D_f) in m
        L_w = 1;                           % Sweep Angle of the Wing (λ_W) in deg
        AR_w = 1;                          % Aspect Ratio (AR_w) in unitless
        twist_w = 1;                       % Wing twist (alpha_t) in deg
        V_C = 1;                           % Cruise Speed (V_C) in m/s
        W_avg = 1;                         % Average weight of aircraft during cruise (W_avg) in kg
        rho = 1.225;                       % Density at cruise (rho) in kg/m^3
        alphaW = 1;                        % Angle of attack of the wing (alphaW) in deg
        alphaF = 1;                        % Angle of attack of the fuselage (alphaF) in deg
        xcg_xw = 1;                        % Distance between centre of gravity and wing 


        % Aerodynamic Coefficients
        CM_af = 1;                         % Wing airfoil section pitching moment coefficient (CM_af) in unitless
        CL_w = 1;                          % Coefficient of lift of the wing (CL_w) in unitless
        CL_alphaW = 1;                     % Wing lift curve slope (CL_αW) in unitless
        Cl_alphaH = 1;                     % Sectional lift coefficient of tailplane airfoil (Cl_Hα) in unitless
        T_ef = 1;                          % Tailplane efficiency factor (Eta, η) in unitless
        horizontalTailVolumeCoefficient = 1; % Horizontal tail volume coefficient (V_barH) in unitless

        % Tail Geometry & Aerodynamics
        sweepAngle_Lambda_ht = 1;          % Sweep angle (Λ_ht) in degrees
        dihedralAngle_Gamma_ht = 1;        % Dihedral angle (Γ_ht) in degrees
        taperRatio_lambdah = 1;            % Taper ratio (λ_h)
        AR_h = 1;                          % Aspect ratio (AR_h)

        % Misc
        KC = 1;                            % Correction Factor (K_C) in unitless 
        
    end 
    
    methods
        % Step 4: Tail arm (L_ht)
        % Use the formula from Equation 6.47
        function val = get.L_ht(obj)
            % Equation: L_ht = K_C * sqrt( (4 * C_bar * S_w * V_barH) / (pi * D_f) )
            val = obj.KC * sqrt(4 * obj.C_bar * obj.S_w * obj.horizontalTailVolumeCoefficient / (pi * obj.Df));
        end
        
        % Step 5: Horizontal tail planform area (S_ht)
        % Apply equation 6.24 to calculate the horizontal tail area.
        function val = get.S_ht(obj)
            % Equation: S_ht = V_barH * C_bar * S_w / L_ht
            val = obj.horizontalTailVolumeCoefficient * obj.C_bar * obj.S_w / obj.L_ht;
        end

        % Step 6: Wing-Fuselage aerodynamic pitching moment coefficient (CM_owf)
        % Calculate the wing-fuselage pitching moment coefficient using equation 6.26.
        function val = get.CM_owf(obj)
            % Equation: CM_owf = CM_af * (AR_w * cos(λ_W)^2) / (AR_w + 2 * cos(λ_W)) + 0.01 * twist_w
            val = obj.CM_af * (obj.AR_w * cos(deg2rad(obj.L_w))^2) / (obj.AR_w + 2 * cos(deg2rad(obj.L_w))) + 0.01 * obj.twist_w;
        end
        
        % Step 7: Cruise lift coefficient (CL_c)
        % Use equation 6.27 to calculate the cruise lift coefficient.
        function val = get.CL_c(obj)
            % Equation: CL_c = (2 * W_avg) / (rho * V_C^2 * S_w)
            val = (2 * obj.W_avg) / (obj.rho * obj.V_C^2 * obj.S_w);
        end
        
        % Step 8: Desired horizontal tail lift coefficient at cruise (CL_h)
        % Use equation 6.29 to calculate the desired lift coefficient for the horizontal tail during cruise via wing properties.
        function val = get.CL_h(obj)
            % Equation: CL_h = (CM_owf + CL_c * L_ht) / (V_barH * T_ef)
            val = (obj.CM_owf + obj.CL_c * obj.L_ht) / (obj.horizontalTailVolumeCoefficient * obj.T_ef);
        end
        
        % Step 12: Horizontal tail lift curve slope (CL_alphaH)
        % Use equation 6.57 to determine the lift curve slope for the horizontal tail.
        function val = get.CL_alphaH(obj)
            % Equation: CL_alphaH = Cl_alphaH / (1 + Cl_alphaH / (pi * AR_h))
            val = obj.Cl_alphaH / (1 + obj.Cl_alphaH / (pi * obj.AR_h));
        end
        
        % Step 13: Horizontal tail incidence angle at cruise (incidence_i_ht)
        % Calculate the angle of attack for the horizontal tail during cruise using equation 6.51.
        function val = get.incidence_i_ht(obj)
            % Equation: incidence_i_ht = CL_h / CL_alphaH  
            val = obj.CL_h / obj.CL_alphaH; 
        end
        
        % Step 14: Downwash angle (epsilon)
        % Calculate the downwash angle at the tail using equation 6.54
        function val = get.epsilon(obj)
            % Equation: epsilon = 2 * CL_w / (pi * AR_w) + 2 * CL_alphaW * alphaW / (pi * AR_w)
            val = 2 * obj.CL_w / (pi * obj.AR_w) + 2 * obj.CL_alphaW * obj.alphaW / (pi * obj.AR_w);
        end 
        
        % Step 15: Effective angle of attack of the tailplane (alphaH_eff)
        % Use equation 6.53 to calculate the incidence angle for the horizontal tail.
        function val = get.alphaH_eff(obj)
            % Equation: alphaH_eff = alphaF + incidence_i_ht - epsilon
            val = obj.alphaF + obj.incidence_i_ht - obj.epsilon;
        end

        % Step 16: Tail dimensions calculations
        % Apply equations 6.63 through 6.66 to calculate the tail span, root chord, tip chord, and mean aerodynamic chord.
        function val = get.TailDimensions(obj)
            % Define symbolic variables for the tail geometry
            syms B_ht C_ht_tip C_ht_root meanAerodynamicChord_ht

            % Equation for AR_h: AR_h = B_ht / meanAerodynamicChord_ht
            AR_eqn = obj.AR_h == B_ht / meanAerodynamicChord_ht;
            
            % Equation for taper ratio: lambda_h = C_ht_tip / C_ht_root
            lambda_eqn = obj.taperRatio_lambdah == C_ht_tip / C_ht_root;
            
            % Equation for mean aerodynamic chord of the horizontal tail:
            % meanAerodynamicChord_ht = (2/3) * C_ht_root * (1 + lambda_h + lambda_h^2) / (1 + lambda_h)
            meanAerodynamicChord_eqn = (2/3) * C_ht_root * (1 + obj.taperRatio_lambdah + obj.taperRatio_lambdah^2) / (1 + obj.taperRatio_lambdah) == meanAerodynamicChord_ht;
            
            % Equation for planform area of the horizontal tail: S_ht = B_ht * meanAerodynamicChord_ht
            S_ht_eqn = B_ht * meanAerodynamicChord_ht == obj.S_ht;
            
            % Solve the system of equations
            sol = solve([meanAerodynamicChord_eqn, S_ht_eqn, AR_eqn, lambda_eqn], [B_ht, C_ht_tip, C_ht_root, meanAerodynamicChord_ht]);
            
            % Store the solutions in a structure
            val = struct();
            val.B_ht = double(sol.B_ht);
            val.C_ht_tip = double(sol.C_ht_tip);
            val.C_ht_root = double(sol.C_ht_root);
            val.meanAerodynamicChord_ht = double(sol.meanAerodynamicChord_ht);
            
            % Display the results for debugging
            disp('B_ht:')
            disp(val.B_ht)
            
            disp('C_ht_tip:')
            disp(val.C_ht_tip)
            
            disp('C_ht_root:')
            disp(val.C_ht_root)
            
            disp('meanAerodynamicChord_ht:')
            disp(val.meanAerodynamicChord_ht)
        end

        % Step 17: Generated lift coefficient at cruise (CLH_c)
        function val = get.CLH_c(~)
            % Placeholder value for the generated horizontal tail lift coefficient at cruise
            val = 2; % Placeholder value
        end
        
        % Testing: Lift and stability checks:
        function resultLift = checkConditionLift(obj)
            % If the generated lift coefficient does not match the required value, adjust the tail incidence angle.
            if obj.CL_h < obj.CLH_c
                resultLift = 'Lift generated by tailplane at cruise is greater than needed for trimmed flight, change incidence angle or re-evaluate other parameters';
            elseif obj.CL_h > obj.CLH_c
                resultLift = 'Lift generated by tailplane at cruise is less than needed for trimmed flight, change incidence angle or re-evaluate other parameters';
            else
                resultLift = 'Lift is exactly to sustain trimmed flight, well done';
            end
            disp(resultLift);
            disp(['Trimmed flight CL required: ', num2str(obj.CL_h)]);
            disp(['Actual CL: ', num2str(obj.CLH_c)]);
            fprintf('\n');
        end

        function resultLongitudinal = checkConditionLongitudinal(obj)
            % Ensure that the Cmα derivative is negative for stabilizing contribution. If not, redesign the tail.
            CMa = obj.CL_alphaW*(obj.xcg_xw/obj.C_bar) - obj.CL_alphaH*obj.T_ef*obj.S_ht/obj.S_w *(1 - 2 * obj.CL_alphaW / (pi * obj.AR_w))*(obj.L_ht/obj.C_bar) ; % Placeholder for actual calculation
            if CMa < 0
                resultLongitudinal = 'Derivative condition for horizontal tailplane is fulfilled';
            else
                resultLongitudinal = 'Derivative condition for horizontal tailplane is NOT fulfilled';
            end
            disp(resultLongitudinal);
            disp(['Actual CM_alpha value: ', num2str(CMa)])
            fprintf('\n');
        end
        
        function resultOverall = checkOverall(obj)
            resultLift = obj.checkConditionLift();
            resultLongitudinal = obj.checkConditionLongitudinal();
            if strcmp(resultLift, 'Lift is sufficient to sustain trimmed flight') && ...
               strcmp(resultLongitudinal, 'Longitudinal stability requirements are fulfilled')
                resultOverall = 'Both derivative and lift conditions are fulfilled';
            elseif ~strcmp(resultLongitudinal, 'Longitudinal stability requirements are fulfilled') && ...
                   strcmp(resultLift, 'Lift is sufficient to sustain trimmed flight')
                resultOverall = 'Derivative condition not fulfilled';
            elseif strcmp(resultLongitudinal, 'Longitudinal stability requirements are fulfilled') && ...
                   ~strcmp(resultLift, 'Lift is sufficient to sustain trimmed flight')
                resultOverall = 'Lift condition not fulfilled';
            else
                resultOverall = 'Both derivative and lift conditions are NOT fulfilled';
            end
            disp(resultOverall);
            fprintf('\n');
        end
        
        % More checks - to be done later / qualitatively

        % Step 19: Check horizontal tail stall
        % Verify that the horizontal tail does not stall under the given conditions 
        % >Done via qualitative empirical check

        % Step 21: Analyze dynamic longitudinal stability
        % Perform an analysis of the dynamic longitudinal stability of the aircraft. If the design requirements are not satisfied, redesign the tail.
        % >Later consideration as it needs whole aricraft consideration
        
        % Step 22: Optimize horizontal tail
        % Optimize the horizontal tail design based on the above calculations and analysis to meet stability and performance requirements.
        % >Qualitative rather than quantitative

        % Note, elevator sizing is done in G17 via vibes hence no calcs are
        % needed necessarily
    end
end
