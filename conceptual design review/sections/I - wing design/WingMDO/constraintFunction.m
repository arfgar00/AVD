function [c, ceq] = constraintFunction(x)
    % Existing constraints
    optIdx = wingMDO(x);
    c = [-optIdx];  % Ensures optIdx >= 0
    
    % Design variables
    cr = 12;
    ck = x(1);
    ct = x(2);
    s = 65/2;
    yk = x(5);
    
    % Inequality constraints (c <= 0)
    c = [
        c;                % Existing constraints
        ck - cr;          % Ensures ck <= cr
        ct - ck;          % Ensures ct <= ck
        yk - s            % Ensures yk <= s
    ];
    
    % Desired wing area
    S_desired = 515.2;  % Replace with your target wing area value
    
    % Calculate the wing area
    wing = WingGeometry();
    wing.cr = cr;
    wing.ck = ck;
    wing.ct = ct;
    wing.s = s;
    wing.Lambdain50 = x(5);
    wing.Lambdaout50 = x(6);
    wing.yk = yk;
    wing = wing.calcSref();
    area = wing.SREF;
    
    % Equality constraint (area == S_desired)
    ceq = area - S_desired;
end
