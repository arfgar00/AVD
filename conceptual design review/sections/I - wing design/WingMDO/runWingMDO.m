run("SetParametersMDO.m")

set(groot, 'defaultFigureUnits', 'inches');
set(groot, 'defaultFigurePosition', [1, 1, 6, 4]); % Width = 6 inches, Height = 4 inches

% Set default font settings for all text in plots
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12); % Font size for axis labels and ticks
set(groot, 'defaultTextFontSize', 12); % Font size for text in the plot
set(groot, 'defaultAxesTitleFontSizeMultiplier', 1.2); % Slightly larger title
set(groot, 'defaultLegendFontSize', 11); % Font size for legend text

% Set default line and marker settings
set(groot, 'defaultLineLineWidth', 1.5); % Default line width
set(groot, 'defaultLineMarkerSize', 8); % Default marker size

% Set default box and grid
set(groot, 'defaultAxesBox', 'on'); % Enable box around plots by default
set(groot, 'defaultAxesGridLineStyle', '-'); % Solid grid lines
set(groot, 'defaultAxesXGrid', 'on'); % Enable grid for X axis
set(groot, 'defaultAxesYGrid', 'on'); % Enable grid for Y axis
% Initial Guess of wing
%x = [cr, ck, ct, Lambdaout50, yk, maxtwist in radians]
%x0 = [14, 12, 4, 30*pi/180, 30*pi/180, 12, deg2rad(1.3)];
x0 = [14, 12, 1.7, 30*pi/180, 10, deg2rad(1.3)];

% Show performance of initial wing
wing0 = x2wing(x0);

% Bounds
lb = [5, 5, 1, 10*pi/180, bodyDiameter*1.2, deg2rad(0)];
ub = [23, 15, 3.5, 35*pi/180, 14, deg2rad(10)];

% Setup MDO
options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...       % Choose an appropriate algorithm
    'OptimalityTolerance', 1e-6, ...         % Set Optimality Tolerance
    'FunctionTolerance', 1e-6, ...           % Set Function Tolerance
    'StepTolerance', 1e-12, ...              % Set Step Tolerance
    'MaxIterations', 1000, ...               % Set Max Iterations
    'MaxFunctionEvaluations', 5000, ...      % Set Max Function Evaluations
    'Display', 'iter');                      % Display iteration output
%options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
[x_opt, fval] = fmincon(@wingMDO, x0, [], [], [], [], lb, ub, @constraintFunction, options);
clc
wing_opt =  x2wing(x_opt);
wing_opt = wing_opt.calcSref();

% Plot two wings, before and after
figure(1)
clf;
wing0.plotWing("b")
hold on
wing_opt.plotWing("r")
hold on
yline(bodyDiameter/2)
hold on
yline(bodyDiameter/2, '--', 'FuseLage', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');
hold on
h1 = plot(NaN, NaN, 'r', 'LineWidth', 2, 'DisplayName', 'Optimized');
hold on
h2 = plot(NaN, NaN, 'b', 'LineWidth', 2, 'DisplayName', 'Initial');
legend([h1, h2], {'Optimized', 'Initial'});
ylim([0 bdes/2]);
print(gcf, 'compareWing', '-dpng', '-r300');

% Show performance of opt wing
disp("wing0: ")
showPerformance(x0,bdes, cruise1, Damping, tolerance,bodyDiameter, l_d, Swet_SrefWing, airfoil)
disp("wingopt: ")
showPerformance(x_opt,bdes, cruise1, Damping, tolerance,bodyDiameter, l_d, Swet_SrefWing, airfoil)

disp("Wf_opt - Wf0/ Wf0")
Wf0 = 0.4143*420739.2;
Wf_opt = 0.400*407954.8;
disp((Wf0-Wf_opt) / Wf_opt)

function showPerformance(x, bdes,cruise1, Damping, tolerance,bodyDiameter, l_d, Swet_SrefWing, airfoil)
    mywing = x2wing(x);
    Wto = wingMDO(x);
    bodyLength = bodyDiameter*l_d;
    Swet_SrefBody = pi*bodyDiameter*bodyLength*1.1/mywing.SREF;    
    cruise1 = cruise1.init(mywing.cbar);
    [CL1, CDi1, Cly1, CDprofile1] = LLESwept(mywing, airfoil, cruise1, Damping, tolerance, bodyDiameter);
    CDF1 = CDFfun(cruise1.M, cruise1.Re, mywing, airfoil.x_cm, airfoil.t_c, l_d, Swet_SrefWing, Swet_SrefBody);
    CDW1 = CDWfun(mywing,Cly1,cruise1.M,airfoil.t_c);
    CDtotal1 = CDF1 + CDW1 + CDi1 + CDprofile1;

    disp(["L/D: ", CL1/ CDtotal1])
    disp(["CLclean: ", CL1])
    disp(["Sref: ", mywing.SREF])
    disp(["W/S: ", Wto/mywing.SREF])
    disp(["Wto: ", Wto])
    disp(["Maximum twist(degrees)", mywing.twist_max * 180/pi])
    disp(["Inner sweep angle(degrees)", mywing.Lambdain50 * 180/pi])
    disp(["Outer sweep angle(degrees)", mywing.Lambdaout50 * 180/pi])
    %disp(mywing)
end
save("wing_opt.mat","wing_opt")