% Initial guess
x0 = [12, 8, 3, 30, 31*pi/180, 25*pi/180, 12];
% Design variables bounds
lb = x0.*0.5;
ub = x0.*2;


wing0 = x02wing(x0);
figure(1)
clf;
axis equal
wing0.plotWing()
hold on
% Optimization options
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

% Run the optimization
[x_opt, fval] = fmincon(@wingMDO, x0, [], [], [], [], lb, ub, @constraintFunction, options);

% Display optimal design variables
disp('Optimal design variables:');
disp(['cr (Root chord length) = ', num2str(x_opt(1)), ' meters']);
disp(['ck (Kink chord length) = ', num2str(x_opt(2)), ' meters']);
disp(['ct (Tip chord length) = ', num2str(x_opt(3)), ' meters']);
disp(['s (Total wingspan) = ', num2str(x_opt(4)), ' meters']);
disp(['Lambdain50 (Inboard sweep) = ', num2str(x_opt(5)*180/pi), ' degrees']);
disp(['Lambdaout50 (Outboard sweep) = ', num2str(x_opt(6)*180/pi), ' degrees']);
disp(['yk (Kink position) = ', num2str(x_opt(7)), ' meters']);

% Display optimal performance index
disp(['Optimal performance index: ', num2str(fval)]);

% Validate the wing area constraint
% Compute wing area for optimal design
cr_opt = x_opt(1);
ck_opt = x_opt(2);
ct_opt = x_opt(3);
s_opt = x_opt(4);
yk_opt = x_opt(7);

wing_opt = x02wing(x_opt);

wing_opt.plotWing
wing_opt.SREF