% Create an instance of the HorizontalTail class
horizontalTail = HorizontalTail();


% Initialize the Airfoil object
airfoil = Airfoil;

% Read polar data from the file
airfoil = airfoil.readPolar('/Users/arfredgarcia/Documents/AVD/AVD GitHub Code/AVD/conceptual design review/sections/I - wing design/WingMDO/xf-sc20614-il-1000000.csv');

% Specify the desired Cl
targetCl = 0.5;  % Replace with your desired Cl value

% Ensure Cl is within the range of the data
if targetCl < min(airfoil.Cl) || targetCl > max(airfoil.Cl)
    error('The specified Cl value is outside the range of the data.');
end

% Interpolate to find alpha
alpha = interp1(airfoil.Cl, airfoil.alpha, targetCl);

% Given angle of attack in radians
queryAlpha = alpha; % Replace with your desired alpha value

% Check if alpha is within range
if queryAlpha < min(airfoil.alpha) || queryAlpha > max(airfoil.alpha)
    error('The specified alpha value is outside the range of the data.');
end

% Interpolate to find Cm for the given alpha
Cm_query = interp1(airfoil.alpha, airfoil.Cm, queryAlpha);

% Display the result
fprintf('The pitching moment coefficient (Cm) for alpha = %.2f degrees is Cm = %.4f.\n', queryAlpha * 180 / pi, Cm_query);

twist_max = 1; % Example value, replace with your value
s = 10;        % Example span, replace with your value

% Define the function
f = wing_opt.twistfun

fprintf('break')

% Compute the average
% Define the function
f = @(y) (twist_max / s) * abs(y);

% Compute the average
average_twist = (1 / s) * integral(f, -s/2, s/2);


% Wing Parameters
horizontalTail.C_bar = wing_opt.cbar;                                                           % Mean aerodynamic chord of the wing (m)
horizontalTail.S_w = wing_opt.SREF;                                                             % Wing reference area (m^2)
horizontalTail.L_w = ( wing_opt.Lambdain50 + wing_opt.Lambdaout50 ) / 2;                        % Sweep angle of the wing (rad)
horizontalTail.AR_w = wing_opt.AR;                                                              % Aspect ratio of the wing
horizontalTail.CM_af = Cm_query;                                                                % Wing airfoil section pitching moment coefficient
horizontalTail.twist_w = average_twist;                                                         % Wing twist (rad)
% horizontalTail.CL_w = targetCl;                                                                 % Coefficient of lift of the wing
horizontalTail.CL_alphaW = 5;                                                                   % Lift curve slope for the wing (1/rad)
horizontalTail.alphaW = alpha;                                                                  % Angle of attack of the wing (rad)

wing = WingGeometry();
wing.cr = wing_opt.cr;
wing.ck = wing_opt.ck;
wing.ct = wing_opt.ct;
wing.s = wing_opt.s;
wing.Lambdain50 = wing_opt.Lambdain50;
wing.Lambdaout50 = wing_opt.Lambdaout50;
wing.yk = wing_opt.yk;
wing.twist_max = 3 * pi()/180;
wing = wing.calcSref()
wing.N = 31;
% figure(3)
% clf;
% axis equal
% wing.plotWing()
wing.c_at_y(10);
wing = wing.createStrips();
% wing.plotStrips();
wing = wing.calcSc() 

c = wing_opt.cbar;
cruise = AirCondition()
cruise.M = 0.8;
cruise.h = convlength(39000, 'ft','m');
cruise = cruise.init(wing_opt.cbar)

AeroWing = LLESwept(wing, airfoil, cruise);                                               % [CL, CDi, CLy]







% Tailplane

% average tailplane needs to be 2 % * MAC lower thickness then that of the
% wing ie wing_opt.cbar (currently 8.5) * 0.02 = 0.17m ?
% ie currently have an sc20712 so should be thinner than this
% should also be symmetric

% Initialize the Airfoil object
tailplaneAirfoil = Airfoil;







% Fuselage Parameters
horizontalTail.Df = 6.38;                                                                       % Fuselage diameter (m), taken from latex report of diameter of fuselage
horizontalTail.alphaF = 0;                                                                      % Angle of attack of the fuselage (rad)


% Tailplane Parameters
horizontalTail.KC = 1.1;                                                                        % Tailplane configuration factor, taken from book
horizontalTail.T_ef = 0.85;                                                                     % Tailplane efficiency factor, taken from book
horizontalTail.Cl_alphaH = 6;                                                                   % Sectional lift coefficient of tailplane airfoil -
horizontalTail.horizontalTailVolumeCoefficient = 1.1;                                           % Horizontal Tail Volume Coefficient (unitless) -
horizontalTail.sweepAngle_Lambda_ht = wing_opt.Lambdain50;                                                       % Sweep angle of tailplane (in degrees or radians), taken from A380
horizontalTail.dihedralAngle_Gamma_ht = 20;                                                     % Dihedral angle of tailplane (in degrees or radians), taken from A380
horizontalTail.taperRatio_lambdah = 0.7;                                                        % Taper ratio of tailplane (unitless)
horizontalTail.AR_h = 7.4;                                                                      % Aspect ratio of the tailplane (unitless), taken from A380 
horizontalTail.CL_w = AeroWing(1);

% Cruise Conditions
horizontalTail.V_C = 250;                                                                       % Cruise speed (m/s)
horizontalTail.W_avg = 300000;                                                                  % Average weight of the aircraft during cruise (N)
horizontalTail.rho = 1.225;                                                                     % Air density at cruise (kg/m^3)

% Misc parameters
horizontalTail.xcg_xw = 3;                                                                      % Distance between CG and wing (in meters)

% Defining the tailplane characteristics such that we can define the lift
% coefficient of the wing
myTailplane = WingGeometry();
myTailplane.cr = horizontalTail.TailDimensions.C_ht_root(horizontalTail.TailDimensions.C_ht_root > 0);                                                                     % Root chord (in meters) 
myTailplane.ct = horizontalTail.TailDimensions.C_ht_tip(horizontalTail.TailDimensions.C_ht_tip > 0);                                                                   % Tip chord (in meters)
myTailplane.ck = 0.5 * (myTailplane.cr + myTailplane.ct);                                        % Chord at mid-span (in meters)
myTailplane.s = horizontalTail.TailDimensions.B_ht(horizontalTail.TailDimensions.B_ht > 0)/2;                                                                            % Semi-wingspan
myTailplane.Lambdain50 = horizontalTail.sweepAngle_Lambda_ht * pi / 180;                                                          % Inboard sweep angle at 50% of the tailplane span (in radians)
myTailplane.Lambdaout50 = horizontalTail.sweepAngle_Lambda_ht * pi / 180;                                                         % Outboard sweep angle at 50% of the tailplane span (in radians)
myTailplane.yk = myTailplane.s/2;                                                                         % Kink position
myTailplane.twist_max = 3 * pi()/180;
myTailplane = myTailplane.calcSref();                                                            % Calculate wing reference area
myTailplane.N = 31;             
% figure(3)
% clf;
% axis equal
myTailplane.plotWing()
%myTailplane.c_at_y(9);
myTailplane = myTailplane.createStrips();
myTailplane.plotStrips();
myTailplane = myTailplane.calcSc

AeroHorizontal = LLESwept(myTailplane, airfoil, cruise);                                               % [CL, CDi, CLy]

horizontalTail.CLH_c = AeroHorizontal(1);

% Check if the horizontal tail functions are fulfilled1
OverallResult = horizontalTail.checkOverall();
Structure = horizontalTail.TailDimensions();
    





      % cr          % Root chord length
      % ck          % Kink chord length
      % ct          % Tip chord length
      % s           % semi-wingspan
      % yk          % Spanwise location of the kink
      % Lambdain50  % Inboard sweep angle at 50% chord
      % Lambdaout50 % Outboard sweep angle at 50% chord
      % N           % Strip theory: number of strips on one wing
      % Sc          % Area of each strip
      % stripy      % Mid point y coordinate of each strip, size = 1x(Nin + Nout)
      % AR          % Aspect ratio
      % Sin         % Area of inner section (2 wing)
      % Sout        % Area of outer section (2 wing)
      % cbar        % Mean aerodynamic chord length
      % b           % Span
      % SREF        % Reference area on one side
      % twist_max   % maximum twist in radians
      % cn          % chord length at strip y points
      % twistfun    % Function of twist at y, a linear twist from 0 to twist_max