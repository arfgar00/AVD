classdef controlSurface
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        y_s            % starting point of the control surface normalized by semi-span
        dy_s           % aileron delta y / semispan
        w_c            % Width of control surface normalized by local chord
        color          % color of control surface on plot
        HLDtype        % name of control surface
        area           % area of the control surface
        LambdaHL       % sweep angle of high lift device
        location       % LE or TE
        cf_c           % Extended chord length/chord length
        dClmax         % Sectional dClmax
    end

    methods
        function obj = calcClInc(obj)
            if obj.HLDtype == "Fowler"
                obj.dClmax = 1.3*obj.cf_c;
            elseif obj.HLDtype == "DoubleSlotted"
                obj.dClmax = 1.6*obj.cf_c;
            elseif obj.HLDtype == "TripleSlotted"
                obj.dClmax = 1.9*obj.cf_c;
            elseif obj.HLDtype == "FixedSlotted"
                obj.dClmax = 0.2;
            elseif obj.HLDtype == "LeadingEdgeFlap"
                obj.dClmax = 0.3;
            elseif obj.HLDtype == "Slat"
                obj.dClmax = 0.4*obj.cf_c;
            end
        end
        function obj = calcS(obj, wing)
            %area is a trapezium. 
            % height is dy_s*s
            h = obj.dy_s*wing.s;
            %c1 and c2 are at y_s and y_s + dy_s
            c1 = (obj.w_c + obj.cf_c - 1) * wing.c_at_y(obj.y_s*wing.s);
            c2 = (obj.w_c + obj.cf_c - 1) * wing.c_at_y((obj.y_s+ obj.dy_s)*wing.s);
            obj.area = (c1 + c2)*h/2;
            if obj.location == "TE"
                obj.LambdaHL = wing.Lambdax_c(1,obj.y_s*wing.s);
            else
                obj.LambdaHL = wing.Lambdax_c(0,obj.y_s*wing.s);
            end
        end
        
        function plotCSTE(obj, wing)
            % The aileron is defined by 3 line segments
            % Vertices: P1, P2, P3, P4
            
            % Points on the trailing edge (TE)
            P1 = [wing.TEx(obj.y_s * wing.s), obj.y_s * wing.s];
            P2 = [wing.TEx((obj.y_s + obj.dy_s) * wing.s), (obj.y_s + obj.dy_s) * wing.s];
            
            % Points on the wing
            P3 = [wing.TEx(obj.y_s * wing.s) - obj.w_c * wing.c_at_y(obj.y_s * wing.s), obj.y_s * wing.s];
            P4 = [wing.TEx((obj.y_s + obj.dy_s) * wing.s) - obj.w_c * wing.c_at_y((obj.y_s + obj.dy_s) * wing.s), (obj.y_s + obj.dy_s) * wing.s];
            
            % Points out of the wing, represents extension of the flap
            P5 = [wing.TEx(obj.y_s * wing.s) + (obj.cf_c - 1) * wing.c_at_y(obj.y_s * wing.s), obj.y_s * wing.s];
            P6 = [wing.TEx((obj.y_s + obj.dy_s) * wing.s) + (obj.cf_c - 1) * wing.c_at_y((obj.y_s + obj.dy_s) * wing.s), (obj.y_s + obj.dy_s) * wing.s];

            % Plotting the aileron
            hold on
            plot([P1(1), P2(1)], [P1(2), P2(2)], obj.color);
            plot([P1(1), P3(1)], [P1(2), P3(2)], obj.color);
            plot([P2(1), P4(1)], [P2(2), P4(2)], obj.color);
            plot([P3(1), P4(1)], [P3(2), P4(2)], obj.color);
            if obj.cf_c > 0
                hold on
                plot([P1(1), P5(1)], [P1(2), P5(2)], obj.color,'LineStyle','--');
                plot([P5(1), P6(1)], [P5(2), P6(2)], obj.color,'LineStyle','--');
                plot([P6(1), P2(1)], [P6(2), P2(2)], obj.color,'LineStyle','--');
            end
            
            % Compute the mean point between P1 and P2
            P_mean = (P1 + P2) / 2;
            
            % Define the end point for the horizontal line
            L_line = 0.2 * wing.cbar; % Adjust the length as needed
            P_end = [P_mean(1) + L_line, P_mean(2)];
            
            % Plot the horizontal line from P_mean towards the positive x-direction
            plot([P_mean(1), P_end(1)], [P_mean(2), P_end(2)],'k--');
            
            % Add the label with obj.HLDtype
            name_latex = strrep(obj.HLDtype, '_', '\_');
            text(P_end(1), P_end(2), ['\textrm{' name_latex '}'], ...
                'Interpreter', 'latex', ...
                'FontSize', 20, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'bottom');
        end

        function plotCSLE(obj, wing)
            % The aileron is defined by 3 line segments
            % Vertices: P1, P2, P3, P4
            
            % Points on the leading edge (LE)
            P1 = [wing.LEx(obj.y_s * wing.s), obj.y_s * wing.s];
            P2 = [wing.LEx((obj.y_s + obj.dy_s) * wing.s), (obj.y_s + obj.dy_s) * wing.s];
            
            % Points on the wing
            P3 = [wing.LEx(obj.y_s * wing.s) + obj.w_c * wing.c_at_y(obj.y_s * wing.s), obj.y_s * wing.s];
            P4 = [wing.LEx((obj.y_s + obj.dy_s) * wing.s) + obj.w_c * wing.c_at_y((obj.y_s + obj.dy_s) * wing.s), (obj.y_s + obj.dy_s) * wing.s];
            
            % Points out of the wing, represents extension of the flap
            P5 = [wing.LEx(obj.y_s * wing.s) - (obj.cf_c - 1) * wing.c_at_y(obj.y_s * wing.s), obj.y_s * wing.s];
            P6 = [wing.LEx((obj.y_s + obj.dy_s) * wing.s) - (obj.cf_c - 1) * wing.c_at_y((obj.y_s + obj.dy_s) * wing.s), (obj.y_s + obj.dy_s) * wing.s];
            
           
            % Plotting the aileron
            hold on
            plot([P1(1), P2(1)], [P1(2), P2(2)], obj.color);
            plot([P1(1), P3(1)], [P1(2), P3(2)], obj.color);
            plot([P2(1), P4(1)], [P2(2), P4(2)], obj.color);
            plot([P3(1), P4(1)], [P3(2), P4(2)], obj.color);
            hold on
            if obj.cf_c > 0
                hold on
                plot([P1(1), P5(1)], [P1(2), P5(2)], obj.color,'LineStyle','--');
                plot([P5(1), P6(1)], [P5(2), P6(2)], obj.color,'LineStyle','--');
                plot([P6(1), P2(1)], [P6(2), P2(2)], obj.color,'LineStyle','--');
            end
            % Compute the mean point between P1 and P2
            P_mean = (P1 + P2) / 2;
            
            % Define the end point for the horizontal line
            L_line = -0.5 * wing.cbar; % Adjust the length as needed
            P_end = [P_mean(1) + L_line, P_mean(2)];
            
            % Plot the horizontal line from P_mean towards the positive x-direction
            plot([P_mean(1), P_end(1)], [P_mean(2), P_end(2)],'k--');
            
            % Add the label with obj.HLDtype
            name_latex = strrep(obj.HLDtype, '_', '\_');
            text(P_end(1), P_end(2), ['\textrm{' name_latex '}'], ...
                'Interpreter', 'latex', ...
                'FontSize', 20, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'bottom');
        end

    end
end