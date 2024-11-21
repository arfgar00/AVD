function [c, ceq] = constraintFunction(x)
    run("SetParametersMDO.m")

    %Wto = wingMDO(x);
    wing = x2wing(x);
    Sref = wing.SREF;

    %W_S = Wto / Sref;
    
    % Inequality constraints (c <= 0)
    c = [
        wing.ck - wing.cr;          % Ensures ck <= cr
        wing.ct - wing.ck;          % Ensures ct <= ck
        Srefdes - Sref      % Sref > Srefdes
    ];
    
    % Equality constraint eg. (area == S_desired)
    %ceq = Sref - S_desired;
    ceq = [];
end
