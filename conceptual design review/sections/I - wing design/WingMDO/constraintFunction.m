function [c, ceq] = constraintFunction(x)
    % Inequality constraints (c <= 0)
    % Existing inequality constraints (if any)
    optIdx = wingMDO(x);
    c = [-optIdx];  % Ensures optIdx >= 0 ];
    
    % Desired wing area
    S_desired = 515.2;  % Replace with your target wing area value
    
    % Calculate the wing area
    wing = WingGeometry();
    wing.cr = x(1);
    wing.ck = x(2);
    wing.ct = x(3);
    wing.s = x(4);
    wing.Lambdain50 = x(5);
    wing.Lambdaout50 = x(6);
    wing.yk = x(7);
    wing = wing.calcSref();
    area = wing.SREF;
    
    % Equality constraint (area == S_desired)
    ceq = area - S_desired;
end