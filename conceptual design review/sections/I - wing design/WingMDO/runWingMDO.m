clear
run("SetParametersMDO.m")

% Initial Guess of wing
%x = [cr, ck, ct, Lambdain50, Lambdaout50, yk, maxtwist in radians]
x0 = [13, 10, 3.5, 30*pi/180, 30*pi/180, 12, deg2rad(1.3)];

% Show performance of initial wing
wing0 = x2wing(x0, bdes);

% Bounds
lb = [13, 10, 3.5, 20*pi/180, 20*pi/180, 8, deg2rad(0)];
ub = [17, 15, 5, 35*pi/180, 35*pi/180, 14, deg2rad(10)];

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
wing_opt =  x2wing(x_opt, bdes);
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
h1 = plot(NaN, NaN, 'r', 'LineWidth', 2, 'DisplayName', 'Optimized');
hold on
h2 = plot(NaN, NaN, 'b', 'LineWidth', 2, 'DisplayName', 'Initial');
legend([h1, h2], {'Optimized', 'Initial'});
ylim([0 bdes/2]);

% Show performance of opt wing
disp("wing0: ")
showPerformance(x0,bdes, cruise1, Damping, tolerance,bodyDiameter, l_d, Swet_SrefWing, airfoil)
disp("wingopt: ")
showPerformance(x_opt,bdes, cruise1, Damping, tolerance,bodyDiameter, l_d, Swet_SrefWing, airfoil)

function showPerformance(x, bdes,cruise1, Damping, tolerance,bodyDiameter, l_d, Swet_SrefWing, airfoil)
    mywing = x2wing(x,bdes);
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