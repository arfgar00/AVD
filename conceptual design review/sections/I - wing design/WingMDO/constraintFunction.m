function [c, ceq] = constraintFunction(x)
    run("SetParametersMDO.m")

    Wto = wingMDO(x);
   
    % Design variables
    cr = x(1);
    ck = x(2);
    ct = x(3);
    yk = x(6);
   
    % Calculate the wing area
    wing = WingGeometry();
    wing.cr = cr;
    wing.ck = ck;
    wing.ct = ct;
    wing.s = bdes/2;
    wing.Lambdain50 = x(4);
    wing.Lambdaout50 = x(5);
    wing.yk = yk;
    wing = wing.calcSref();
    Sref = wing.SREF;

    %W_S = Wto / Sref;
    
    % Inequality constraints (c <= 0)
    c = [
        -Wto              % Weight must be positive
        ck - cr;          % Ensures ck <= cr
        ct - ck;          % Ensures ct <= ck
        Srefdes - Sref      % Sref > Srefdes
    ];
    
    % Equality constraint eg. (area == S_desired)
    %ceq = Sref - S_desired;
    ceq = [];
end
