function [T0_W0] = Constraint_Cruise(h, V_M, W0_Sref, W_W0, AR, e, CD_0)
% This function is used to calculate the Empty Fuel
% Fraction of the Wide Body Commercial Airliner from
% given mission profile
%
% Wf_W0(h_cruise1, h_cruise2, LD_max, C)
%
% outputs: T0_W0   - sea level thrust to weight ratio
% inputs:  h       - cruise altitude (ft)
%          V_M     - cruise velocity (Mach)
%          W0_Sref - guessed Wing Loading (kg)
%          W_W0    - weight fraction at the beginning of segment
%          AR      - aspect ratio
%          e       - Oswald/span efficiency factor
%          CD_0    - Drag ceofficient at zero lift

g = 9.80665; % Standard acceleration of gravity (m/s^2)

T0_W0 = @(alpha, beta, dh_dt, dV_dt, W0_Sref, AR, e, V, CD0, n, rho) ...
    alpha ./ beta .* (1 ./ V .* dh_dt + 1 ./ g .* dV_dt + ...
    (0.5 .* rho .* V .^ 2 .* CD0) ./ (alpha .* W0_Sref) + ...
    (alpha .* n .^ 2 .* W0_Sref) ./ ...
    (0.5 .* rho .* V .^ 2 .* pi .* AR .* e));

% Basic parameters
h = convlength(h, 'ft','m');

[T, c, P, rho] = atmosisa(h);

V = V_M .* c;  % Cruise velocity (m/s)

beta = T_T0(h);

T0_W0 = T0_W0(W_W0, beta, 0.0, 0.0, W0_Sref, AR, e, V, CD_0, 1, rho);

end
