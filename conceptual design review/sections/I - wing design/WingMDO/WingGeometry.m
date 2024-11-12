classdef WingGeometry
   properties
      cr          % Root chord length
      ck          % Kink chord length
      ct          % Tip chord length
      s           % semi-wingspan
      yk          % Spanwise location of the kink
      Lambdain50  % Inboard sweep angle at 50% chord
      Lambdaout50 % Outboard sweep angle at 50% chord
      N           % Strip theory: number of strips on one wing
      Sc          % Area of each strip
      stripy      % Mid point y coordinate of each strip, size = 1x(Nin + Nout)
      AR          % Aspect ratio
      Sin         % Area of inner section (2 wing)
      Sout        % Area of outer section (2 wing)
      cbar        % Mean aerodynamic chord length
      b           % Span
      SREF        % Reference area on one side
      twist_max   % maximum twist in radians
      cn          % chord length at strip y points
      twistfun    % Function of twist at y, a linear twist from 0 to twist_max
   end
   properties
        dY          % dy of each strip
   end
   
   methods
        % Method to calculate the reference area and aspect ratio
        function obj = calcSref(obj)
             obj.b = 2*obj.s;
             obj.SREF = ((obj.cr + obj.ck) * (obj.yk) / 2 + (obj.ct + obj.ck) * (obj.s - obj.yk) / 2)*2;
             obj.AR = (obj.b)^2 / obj.SREF;
             obj.Sin = 2*(obj.cr + obj.ck)*obj.yk/2;
             obj.Sout = 2*(obj.ck + obj.ct)*(obj.s - obj.yk)/2;
             obj.cbar = ((obj.ck + obj.cr)/2 + (obj.ck + obj.ct)/2)/2;
             obj.twistfun = @(y) obj.twist_max./obj.s.*abs(y);
        end
        % Get sweep angle at coordinate (x_c,y)
         function Lambda = Lambdax_c(obj,x_c,y)
             y = abs(y);
            if y >= 0 && y <= obj.yk
                dy = obj.yk;
                x2 = 0.5*obj.cr + dy*tan(obj.Lambdain50) - (0.5 - x_c)*obj.ck;
                x1 = x_c*obj.cr;
                Lambda = atan((x2 - x1)/dy);
            elseif y > obj.yk && y <= obj.s
                dy = obj.s - obj.yk;
                x2 = 0.5*obj.ck + dy*tan(obj.Lambdaout50) - (0.5 - x_c)*obj.ct;
                x1 = x_c*obj.ck;
                Lambda = atan((x2 - x1)/dy);
            else
                error('Lambdax_c:InvalidY', 'Input y = %f is outside the wing semi-span [0, %f]', y, obj.s);
            end
        end
        % Get chord length at y
        function c = c_at_y(obj,y)
            y = abs(y);
            %sweep angle at x/c = 0, for first section
            Lambdain0 = obj.Lambdax_c(0,y);
            %sweep angle at x/c = 1, for first section
            Lambdain1 = obj.Lambdax_c(1,y);
            %sweep angle at x/c = 0, for second section
            Lambdaout0 = obj.Lambdax_c(0,y);
            %sweep angle at x/c = 1, for second section
            Lambdaout1 = obj.Lambdax_c(1,y);
            if y <= obj.yk
                c = obj.cr + y*(tan(Lambdain1) - tan(Lambdain0));
            elseif y > obj.yk && y <= obj.s
                c = obj.ck + (y-obj.yk)*(tan(Lambdaout1) - tan(Lambdaout0));
            end
        end
        % Plot Wing, with x_cm shown
        function plotWing(obj)
            % Plot the original wing (top half)
            plot([0 obj.cr], [0 0], "b")
            hold on
        
            xk0 = obj.yk * tan(obj.Lambdax_c(0, 0));
            xk1 = obj.cr + obj.yk * tan(obj.Lambdax_c(1, 0));
        
            xt0 = xk0 + (obj.s - obj.yk) * tan(obj.Lambdax_c(0, obj.s));
            xt1 = xk1 + (obj.s - obj.yk) * tan(obj.Lambdax_c(1, obj.s));
        
            plot([0 xk0], [0 obj.yk], "b")
            plot([obj.cr xk1], [0 obj.yk], "b")
            plot([xk0 xk1], [obj.yk obj.yk], "b")
            
            plot([xk0 xt0], [obj.yk obj.s], "b")
            plot([xk1 xt1], [obj.yk obj.s], "b")
            plot([xt0 xt1], [obj.s obj.s], "b")
        
            % Plot mid-line sweep angle for top half
            xk5 = 0.5 * obj.cr + obj.yk * tan(obj.Lambdax_c(0.5, 0));
            xt5 = xk5 + (obj.s - obj.yk) * tan(obj.Lambdax_c(0.5, obj.s));
            plot([0.5 * obj.cr xk5], [0 obj.yk], "--b")
            plot([xk5 xt5], [obj.yk obj.s], "--b")
        
            % Mirror the wing about the x-axis (bottom half)
            plot([0 obj.cr], [0 0], "b")
        
            plot([0 xk0], [0 -obj.yk], "b")
            plot([obj.cr xk1], [0 -obj.yk], "b")
            plot([xk0 xk1], [-obj.yk -obj.yk], "b")
            
            plot([xk0 xt0], [-obj.yk -obj.s], "b")
            plot([xk1 xt1], [-obj.yk -obj.s], "b")
            plot([xt0 xt1], [-obj.s -obj.s], "b")
        
            % Plot mid-line sweep angle for bottom half
            plot([0.5 * obj.cr xk5], [0 -obj.yk], "--b")
            plot([xk5 xt5], [-obj.yk -obj.s], "--b")
        
            % Finalize the plot
            xlabel('X-axis');
            ylabel('Y-axis');
            title('Full Wing Plot (Mirrored about X-axis)');
            axis equal;
            hold off;
        end

        % Create strips on the wing
        function obj =  createStrips(obj)
            if rem(obj.N,2) == 0
                error("N must be an odd number")
            else
                obj.stripy = linspace(-obj.s,obj.s,obj.N);
                obj.dY = mean(diff(obj.stripy))/2;
                for i = 1:length(obj.stripy)
                    obj.cn(i) = obj.c_at_y(abs(obj.stripy(i)));
                end
            end
        end
        % Given y, find coordinate of x on leading edge
        function x = LEx(obj,y)
            y = abs(y);
            xk0 = obj.yk*tan(obj.Lambdax_c(0,0));
            if y <= obj.yk
                x = y*tan(obj.Lambdax_c(0,0));
            elseif y > obj.yk && y <= obj.s
                x = xk0 + (y - obj.yk)*tan(obj.Lambdax_c(0,obj.s));
            end
        end
        % Plot strips on the wing
        function plotStrips(obj)
            for i = 1:length(obj.stripy)
                y = obj.stripy(i);
                plot([obj.LEx(y) obj.LEx(y) + obj.c_at_y(y)], [y y], "--m")
                hold on
                
                y_minus_dY = max(-obj.s, y - obj.dY); % Ensure y within semi span
                plot([obj.LEx(y_minus_dY) obj.LEx(y_minus_dY) + obj.c_at_y(y_minus_dY)], [y_minus_dY, y_minus_dY], "m")
                hold on
                
                y_plus_dY = min(obj.s, y + obj.dY); % Ensure y within semi span
                plot([obj.LEx(y_plus_dY) obj.LEx(y_plus_dY) + obj.c_at_y(y_plus_dY)], [y_plus_dY, y_plus_dY], "m")
                hold on
            end
        end
        %Calculate area of each strip
        function obj = calcSc(obj)
            for i = 1:length(obj.stripy)
                c1 = obj.c_at_y(min(obj.stripy(i) + obj.dY, obj.s));
                c2 = obj.c_at_y(max(-obj.s, obj.stripy(i)-obj.dY));
                
                %Handle boundary
                if i == 1
                    obj.Sc(i) = (c1 + obj.ct)*obj.dY/2;
                elseif i == length(obj.stripy)
                    obj.Sc(i) = (obj.ct + c2)*obj.dY/2;
                elseif  c1 > obj.yk && c2 < obj.yk %Handle kink
                    A1 = (c1 + obj.yk)*(c1 - obj.yk)/2;
                    A2 = (c2 + obj.yk)*(obj.yk - c2)/2;
                    obj.Sc(i) = A1 + A2;
                else
                    obj.Sc(i) = (c1 + c2)*2*obj.dY/2;
                end
            end
        end
        % function of x coordinate at 0.25c, f(y)
        function x = f(obj,yArray)
            yArray = abs(yArray);
            for i = 1:length(yArray)
                y = yArray(i);
                Lambdain = obj.Lambdax_c(0.25,0);
                Lambdaout = obj.Lambdax_c(0.25,obj.s);
                xk = 0.25*obj.cr + obj.yk*tan(Lambdain);
                if y < obj.yk
                    x(i) = 0.25*obj.cr + y * tan(Lambdain);
                elseif (y >= obj.yk) && (y <= obj.s)
                    x(i) = xk + (y - obj.yk) * tan(Lambdaout);
                else
                   error('f:Invalid_Y', 'Input y = %f is outside the wing semi-span [0, %f]', y, obj.s);
                end
            end
       end
   end
end