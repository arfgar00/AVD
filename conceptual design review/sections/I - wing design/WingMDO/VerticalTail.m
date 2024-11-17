classdef VerticalTail
    properties (Dependent)
        % Geometric and Position Properties
        VTailDimensions                   % Structure containing bv Ct_v Cr_v MACv
        Sv                                % Planform area of the vertical tail (square meters)
        lvt                               % Tail arm, distance from CG to vertical tail aerodynamic center (meters)
        % Aerodynamic and Shape Properties
        ARv                               % Aspect ratio of the vertical tail
        lambda_v                          % Taper ratio (λ_v) of the vertical tail
        sweepAngle_v                      % Sweep angle (Λ_v) of the vertical tail (degrees)
        dihedralAngle_v                   % Dihedral angle (Γ_v) of the vertical tail (degrees)
        incidence_v                       % Incidence angle (iv) of the vertical tail (degrees)
    end

    properties (SetAccess = public)
        C_bar = 1;                        % Wing Mean Aerodynamic Chord (C_bar) of the vertical tail (meters)
        S_w = 1;                          % Wing Reference Area (S) of the vertical tail (square meters)
        VerticalTailVolumeCoefficient=1;  % Vertical tail volume coefficient (V_barV) in unitless 
        lv = 1;                           % Tail moment arm (l_v)
        % airfoilSection                  % Airfoil section used for the vertical tail
        % location                        % Vertical tail location (relative to aircraft origin)
    end

    methods
        % Vertical Tail Design Steps
        
        % % Step 2: Select the vertical tail volume coefficient, Vv
        % function VerticalTailVolumeCoefficient = get.VerticalTailVolumeCoefficient(~)
        %     VerticalTailVolumeCoefficient = 0.9;
        % end
        % 
        % % Step 3: Assume the vertical tail moment arm (lv) as equal to the horizontal tail moment arm (l)
        % function lv = get.lv(~)
        %     lv = 0.9;
        % end        
        
        % Step 4: Calculate vertical tail planform area, Sv
        function Sv = get.Sv(~)
            Sv = obj.VerticalTailVolumeCoefficient * obj.S_w * obj.bv / obj.lv;
        end
        
        % Step 6: Select vertical tail aspect ratio, ARv
        function ARv = get.ARv(~)
            ARv = obj.bv / obj.MACv;
        end
        
        % Step 7: Select vertical tail taper ratio, λv
        function lambda_v = get.lambda_v(~)
            lambda_v = obj.Ct_v / obj.Cr_v;
        end
        
        % % Step 8: Determine the vertical tail incidence angle
        % function incidence_v = get.incidence_v(~)
        %     incidence_v = 0.1;
        % end
        % 
        % % Step 9: Determine the vertical tail sweep angle
        % function sweepAngle_v = get.sweepAngle_v(~)
        %     sweepAngle_v = 0.1;
        % end
        % 
        % % Step 10: Determine the vertical tail dihedral angle
        % function dihedralAngle_v = get.dihedralAngle_v(~)
        %     dihedralAngle_v = 0.1;
        % end

        % Step 11: Calculate vertical tail span (bv), root chord (Cr_v), tip chord (Ct_v),
        % and mean aerodynamic chord (MACv)
        function VTailDimensions = get.VTailDimensions(~)
            % Define symbolic variables
            syms bv Ct_v Cr_v MACv

            % Define the aspect ratio as a symbolic relation
            AR_eqn = obj.ARv == bv / MACv;  % Adjusted to use obj.ARv

            % Define the taper ratio as a symbolic relation
            lambda_eqn = obj.lambda_v == Ct_v / Cr_v;  % Adjusted to use obj.lambda_v

            % Equation for MACv (mean aerodynamic chord of the vertical tail)
            meanAerodynamicChord_eqn = (2/3) * Cr_v * (1 + obj.lambda_v + obj.lambda_v^2) / (1 + obj.lambda_v) == MACv;

            % Equation for S_v (planform area of the vertical tail)
            S_v_eqn = bv * MACv == obj.Sv;

            % Solve the system of equations for bv, Ct_v, Cr_v, MACv
            sol = solve([meanAerodynamicChord_eqn, S_v_eqn, AR_eqn, lambda_eqn], [bv, Ct_v, Cr_v, MACv]);

            % Create a structure to store the solutions
            VTailDimensions = struct();
            VTailDimensions.bv = double(sol.bv);  % Solution for bv (span of the vertical tail)
            VTailDimensions.Ct_v = double(sol.Ct_v);  % Solution for Ct_v (tip chord of the vertical tail)
            VTailDimensions.Cr_v = double(sol.Cr_v);  % Solution for Cr_v (root chord of the vertical tail)
            VTailDimensions.MACv = double(sol.MACv);  % Solution for MACv (mean aerodynamic chord of the vertical tail)

            % Display the structure to check the values
            disp('VTailDimensions structure:');
            disp(VTailDimensions);

            % Display the results for debugging
            disp('B_hv:')
            disp(VTailDimensions.B_hv)
            
            disp('C_vt_tip:')
            disp(VTailDimensions.Ct_v)
            
            disp('C_vr_root:')
            disp(VTailDimensions.C_rv)
            
            disp('meanAerodynamicChord_hv:')
            disp(VTailDimensions.MACv)
        end

        % Testing the tailplane:   
        
        % Step 15: Analyze directional stability
        % Reference: Section 6.8.1
        % Bunch of derivatives we do not have


        % More checks - to be done later / qualitatively

        % Step 12: Check the spin recovery
        % Ensure design is capable of safe recovery from spins if applicable
        % >This is a qualitative check

        % Step 13: Adjust the location of the vertical tail relative to the horizontal tail by changing lv,
        % to satisfy the spin recovery requirements
        % Reference: Section 6.8.2-2
        % >This is a qualitative check

        % Step 14: Analyze directional trim
        % Reference: Section 6.8.1
        % >This is enforced by having a symmetric airfoil 

        % Step 16: Modify to meet the design requirements
        % Ensure the vertical tail meets all aerodynamic, structural, and stability requirements
        % >This is a qualitative check
        
        % Step 17: Optimize the tail design
        % Further optimize geometry, configuration, and materials for performance
        % Step 16 and Step 17 are more just final checks and does not
        % require calculation
        % >This is a qualitative check

        % Note, rudder sizing is done in G17 via vibes hence no calcs are
        % needed necessarily
    end
end
