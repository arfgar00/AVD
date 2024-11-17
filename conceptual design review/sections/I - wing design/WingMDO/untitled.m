clear
clc

% Create an instance of the HorizontalTail class
horizontalTail = HorizontalTail();

% Optionally, set any properties as needed for the specific conditions you want to test
horizontalTail.KC = 1.1;            % Example of setting properties
horizontalTail.C_bar = 1.5;          % Set mean aerodynamic chord of the wing
horizontalTail.S_w = 150;            % Set wing reference area
horizontalTail.Df = 5;               % Set fuselage diameter
horizontalTail.L_w = 25;             % Set sweep angle of the wing in rad
horizontalTail.AR_w = 8;             % Set aspect ratio of the wing
horizontalTail.CM_af = -0.05;        % Set wing airfoil section pitching moment coefficient
horizontalTail.twist_w = 2;          % Set wing twist in rad
horizontalTail.V_C = 250;            % Set cruise speed in m/s
horizontalTail.W_avg = 180000;       % Set average weight of the aircraft during cruise in N
horizontalTail.rho = 1.225;          % Set density at cruise in kg/m^3
horizontalTail.T_ef = 0.95;          % Set tailplane efficiency factor
horizontalTail.Cl_alphaH = 5.5;      % Set sectional lift coefficient of tailplane airfoil
horizontalTail.CL_w = 0.5;           % Set coefficient of lift of the wing
horizontalTail.CL_alphaW = 4.5;      % Set horizontal lift curve slope for the wing
horizontalTail.alphaW = 5;           % Set angle of attack of the wing in rad
horizontalTail.alphaF = 2;           % Set angle of attack of the fuselage in rad

    epsilonValue = horizontalTail.epsilon;
disp(epsilonValue)




% Check if the horizontal tail functions are fulfilled
OverallResult = horizontalTail.checkOverall();

Structure = horizontalTail.TailDimensions();
    



