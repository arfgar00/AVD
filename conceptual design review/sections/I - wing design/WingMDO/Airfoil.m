classdef Airfoil
    properties
        alpha          % angle of attack original, an array
        Cl             % Cl polar original, an array
        Cd             % Cd polar original, an array
        uppershape        % coordinate of shape, original
        lowershape        % coordinate of shape, original
        upperShapefun  % polyfit of upper surface shape
        lowerShapefun  % polyfit of lower surface shape
        x_c               % linspace(0,1,500), chordwise location
        t_cfun            % thickness/chord function from interpolation
        t_c               % maximum thickness/chord
        x_cm              % x_c location at maximum thickness
        a0               % CL = a0*alpha + b
        b                % CL = a0*alpha + b
        alphaSPos        % positive stall angle of attack
        alphaSNeg        % negative stall angle of attack
    end
    
    methods
        function obj = readPolar(obj,filename)
	        tab = readtable(filename);
            obj.Cl = tab.Cl;
            obj.alpha = tab.Alpha.*pi/180;
            obj.Cd = tab.Cd;
             % Define the range in radians
            alpha_min_deg = -7;  % Minimum angle in degrees
            alpha_max_deg = 7;   % Maximum angle in degrees
            alpha_min = alpha_min_deg * pi / 180;  % Convert to radians
            alpha_max = alpha_max_deg * pi / 180;  % Convert to radians
            
            % Find indices within the specified range
            idx_linear = (obj.alpha >= alpha_min) & (obj.alpha <= alpha_max);
            
            % Extract the linear region data
            alpha_linear = obj.alpha(idx_linear);
            Cl_linear = obj.Cl(idx_linear);
            
            % -------------------------------
            % Step 4: Check for Sufficient Data Points
            % -------------------------------
            
            % Ensure there are enough points for a reliable fit
            if length(alpha_linear) < 2
                error('Not enough data points within the range of -5 to +5 degrees for linear fit.');
            end
            
            % -------------------------------
            % Step 5: Perform Linear Regression (Polyfit)
            % -------------------------------
            
            % Perform linear fit: Cl = a0 * alpha + b
            p = polyfit(alpha_linear, Cl_linear, 1);
            obj.a0 = p(1);  % Lift curve slope (Cl per radian)
            obj.b = p(2);   % Intercept (Cl at alpha = 0)
            [~, maxIndex] = max(obj.Cl);
            obj.alphaSPos = obj.alpha(maxIndex);
            [~, minIndex] = min(obj.Cl);
            obj.alphaSNeg = obj.alpha(minIndex);
        end

        function plotPolar(obj)
            subplot(2,2,1)
            plot(obj.alpha.*180/pi, obj.Cl)
            hold on
            plot(obj.alpha.*180/pi, obj.alpha*obj.a0 + obj.b,'--')
            xlabel("alpha")
            ylabel("Cl")
            subplot(2,2,2)
            plot(obj.alpha.*180/pi, obj.Cd)
            xlabel("alpha")
            ylabel("Cd")
            subplot(2,2,3)
            plot(obj.alpha.*180/pi, obj.Cl./obj.Cd)
            xlabel("alpha")
            ylabel("Cl/Cd")
        end

        function [Cl, Cd] = interpPolar(obj,alpha)
            Cl = interp1(obj.alpha, obj.Cl, alpha);
            Cd = interp1(obj.alpha, obj.Cd, alpha);
        end

        function obj = readShape(obj,filename)
            obj.x_c = linspace(0,1,500);
            % readAirfoilCoordinates Reads airfoil coordinates and splits them into upper and lower surfaces.
            %
            % Syntax:
            %   [upperSurface, lowerSurface] = readAirfoilCoordinates(filename)
            %
            % Inputs:
            %   filename - String, the name of the airfoil data file.
            %
            % Outputs:
            %   upperSurface - Nx2 array containing [x, y] coordinates of the upper surface.
            %   lowerSurface - Mx2 array containing [x, y] coordinates of the lower surface.
            %
            % Example:
            %   [upper, lower] = readAirfoilCoordinates('airfoil_data.txt');
            %   plot(upper(:,1), upper(:,2), 'b-', lower(:,1), lower(:,2), 'r-');
            
                % Open the file
                fid = fopen(filename, 'r');
                if fid == -1
                    error('Cannot open the file: %s', filename);
                end
            
                % Skip the header lines
                headerLine1 = fgetl(fid); % Airfoil name
                headerLine2 = fgetl(fid); % Numbers (e.g., 103. 103.)
            
                % Initialize arrays to hold coordinates
                coordinates = [];
            
                % Read the data until end of file
                while ~feof(fid)
                    line = fgetl(fid);
                    % Skip empty lines
                    if isempty(line)
                        continue;
                    end
                    % Parse the line
                    data = sscanf(line, '%f %f');
                    if length(data) == 2
                        coordinates = [coordinates; data'];
                    end
                end
            
                % Close the file
                fclose(fid);
            
                % Split the coordinates into upper and lower surfaces
                % Find the index where the data repeats (start of lower surface)
                % Assume that the first point is repeated (0, 0) for closed airfoils
                % and the data is ordered from leading edge to trailing edge for upper
                % surface, then trailing edge to leading edge for lower surface.
            
                % Find the index where x repeats from 1.0 back to 0.0
                x_diff = diff(coordinates(:,1));
                idx_reset = find(x_diff < 0, 1, 'first');
            
                if isempty(idx_reset)
                    error('Could not find the splitting point between upper and lower surfaces.');
                end
            
                % Split the coordinates
                upperSurface = coordinates(1:idx_reset, :);
                lowerSurface = coordinates(idx_reset+1:end, :);
            
                % Optional: Ensure that the surfaces are correctly ordered from leading to trailing edge
                % Reverse the lower surface if necessary
                if lowerSurface(1,1) > lowerSurface(end,1)
                    lowerSurface = flipud(lowerSurface);
                end
                obj.uppershape = upperSurface;
                obj.lowershape = lowerSurface;
        end

        function obj = interpShape(obj, n)
            % Fit polynomials to the upper and lower shapes
            p_upper = polyfit(obj.uppershape(:,1), obj.uppershape(:,2), n);
            p_lower = polyfit(obj.lowershape(:,1), obj.lowershape(:,2), n);
            
            % Define the upper and lower shape functions
            obj.upperShapefun = @(x_c) polyval(p_upper, x_c);
            obj.lowerShapefun = @(x_c) polyval(p_lower, x_c);
            
            % Calculate thickness at each x_c point
            thickness = obj.upperShapefun(obj.x_c) - obj.lowerShapefun(obj.x_c);
            
            % Fit a polynomial to the thickness data
            p_thickness = polyfit(obj.x_c, thickness, n);
            
            % Define the thickness function
            obj.t_cfun = @(x_c) polyval(p_thickness, x_c);
            
            % Find the maximum thickness and its location
            thickness_values = obj.t_cfun(obj.x_c);
            [obj.t_c, max_index] = max(thickness_values);
            obj.x_cm = obj.x_c(max_index);
        end
                
        function plotShape(obj)
            scatter(obj.uppershape(:,1), obj.uppershape(:,2),".")
            hold on
            scatter(obj.lowershape(:,1), obj.lowershape(:,2),".")
            hold on
            plot(obj.x_c, obj.upperShapefun(obj.x_c),"--")
            hold on
            plot(obj.x_c, obj.lowerShapefun(obj.x_c),"--")
            axis equal
        end
    end
end