function W_f_0 = ...
    Wf_W0(cruise1, cruise2, LD_cruise1, LD_cruise2, C_cruise)
% This function is used to calculate the Empty Fuel
% Fraction of the Wide Body Commercial Airliner from
% given mission profile
%
% Wf_W0(h_cruise1, h_cruise2, LD_max, C)
%
% outputs: W_f_0     - empty fuel fraction
%          W_i       - fuel fraction at each segment
% inputs:  h_cruise1 - first cruise average altitude (ft)
%          h_cruise1 - second cruise altitude (ft)
%          V_cruise1 - first cruise velocity (Mach)
%          V_cruise2 - second cruise velocity (Mach)
%          LD_cruise1   - L/D for first cruise
%          LD_cruise2   - L/D_for second cruise
%          C_cruise  - cruise SFC (1/s)

h_cruise1 = cruise1.h;
V_cruise1 = cruise1.V;
h_cruise2 = cruise2.h;
V_cruise2 = cruise2.V;

% Empirical Airplane Data
LD_cruise = LD_cruise1;        % L/D_cruise % Raymer p.41
LD_alternate = LD_cruise2; % L/D_alternate % Roskam p.57

% Empirical Engine Data
C_loiter = 0.8 .* C_cruise;     % SFC (1/s) % Raymer p.36
C_alternate = 1.15 .* C_cruise; % SFC (1/s) % 115% cruise SFC % Ref12

% Empirical Wf_W0 data
W_TO = 0.970;   % Warmup and Takeoff
W_C = 0.985; % Climb and Accelerate to cruise
W_D = 0.990; % Descent to land
W_L = 0.992; % Landing and Taxi

% 0-1 % Taxi & Takeoff
W_10 = W_TO;

% 1-2 % Climb and accelerate to cruise
W_21 = W_C;

% 2-3 % Cruise for an air distance of 7500 nmi at Mach 0.83
R_32 = convlength(7500, 'naut mi','m'); % Range (m)
h_32 = convlength(h_cruise1, 'ft','m'); % Cruise altitude (m)
[T_32, c_32, P_32, rho_32] = atmosisa(h_32);
V_32 = V_cruise1 * c_32;                % Velocity (m/s)
LD_32 = LD_cruise;                      % Optimal cruise L/D ratio
W_32 = exp((-R_32 .* C_cruise) ./ (V_32 .* LD_32));

% 3-4 % Descent to land
W_43 = W_D;

% 4-5 Missed approach,
% followed by clmb and accelerate to alternate
W_54 = 1 - (h_cruise2 / h_cruise1) * (1 - W_C);

% 5-6 Cruise to alternate destination 370 km away
%
% Assume cruise indicate velocity 280 knots [5]
%
R_65 = convlength(370, 'km', 'm');  % Range (m)
h_65 = convlength(18000, 'ft','m'); % Alternate altitude (m)
[T_65, c_65, P_65, rho_65] = atmosisa(h_65);
V_65 = V_cruise2 .* c_65;           % Velocity (m/s)
LD_65 = LD_alternate;               % Estimate alternate L/D ratio
W_65 = exp((-R_65 .* C_alternate) ./ (V_65 .* LD_65));

% 6-7 % Loiter at 5000 ft for 45 minutes
E_76 = 45 * 60; % Endurance (s)
LD_76 = LD_cruise2/0.866; % Optimal loiter L/D ratio
W_76 = exp((-E_76 .* C_loiter) ./ LD_76);

% 7-8 % Descent to land
W_87 = 1 - (h_cruise2 / h_cruise1) *  (1 - W_D);

% 8-9 % Landing & Taxi
W_98 = W_L;

% Culmulative fuel fraction
W_f_0 = 1.01 .* ...
(1 - W_10 .* W_21 .* W_32 .* W_43 .* W_54 .* W_65...
.* W_76 .* W_87 .* W_98);

% Store Value into Matrix
W_i = [W_10 W_21 W_32 W_43 W_54 W_65 W_76 W_87 W_98];

end
