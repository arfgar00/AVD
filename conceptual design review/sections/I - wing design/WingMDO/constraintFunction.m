function [c, ceq] = constraintFunction(x)
    global bdes Srefdes
    % Existing constraints
    optIdx = wingMDO(x);
    c = [-optIdx];  % Ensures optIdx >= 0
    
    % Design variables
    cr = x(1);
    ck = x(2);
    ct = x(3);
    s = bdes/2;
    yk = x(6);
   
    % Desired wing area
    S_desired = Srefdes;  % Replace with your target wing area value
    
    % Calculate the wing area
    wing = WingGeometry();
    wing.cr = cr;
    wing.ck = ck;
    wing.ct = ct;
    wing.s = s;
    wing.Lambdain50 = x(4);
    wing.Lambdaout50 = x(5);
    wing.yk = yk;
    wing = wing.calcSref();
    area = wing.SREF;

    % Inequality constraints (c <= 0)
    c = [
        c;                % Existing constraints
        ck - cr;          % Ensures ck <= cr
        ct - ck;          % Ensures ct <= ck
        yk - s;            % Ensures yk <= s
    ];
    
    
    %ceq = [];
    % Equality constraint (area == S_desired)
    ceq = area - S_desired;
end
