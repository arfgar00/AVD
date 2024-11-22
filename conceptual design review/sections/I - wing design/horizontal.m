clc

Vh = 1.1
C_bar = 7.5418
Sref = 351.0294 
L_ht = 33.4
Sh = Vh * C_bar * Sref / L_ht

b = 65
Vv = 0.09;
L_vt = 32.4;
Sv = Vv * b * Sref / L_vt

% Define variables
AR_h = 8;                  % Aspect Ratio of the horizontal tail
taperRatio_lambdah = 0.5;    % Taper ratio of the horizontal tail
S_ht = Sh;                   % Planform area of the horizontal tail in m^2

% Solve for the mean aerodynamic chord
syms meanAerodynamicChord_ht C_ht_root B_ht C_ht_tip

% Equation for the mean aerodynamic chord
meanAerodynamicChord_eqn = (2/3) * C_ht_root * ...
    (1 + taperRatio_lambdah + taperRatio_lambdah^2) / ...
    (1 + taperRatio_lambdah) == meanAerodynamicChord_ht;

% Equation for the span of the horizontal tail
AR_eqn = AR_h == B_ht / meanAerodynamicChord_ht;

% Equation for the taper ratio
lambda_eqn = taperRatio_lambdah == C_ht_tip / C_ht_root;

% Equation for the planform area
S_ht_eqn = B_ht * meanAerodynamicChord_ht == S_ht;

% Solve the system of equations
sol = solve([meanAerodynamicChord_eqn, S_ht_eqn, AR_eqn, lambda_eqn], ...
    [B_ht, C_ht_tip, C_ht_root, meanAerodynamicChord_ht]);

% Extract solutions and filter for positive values
B_ht = double(sol.B_ht);
C_ht_tip = double(sol.C_ht_tip);
C_ht_root = double(sol.C_ht_root);
meanAerodynamicChord_ht = double(sol.meanAerodynamicChord_ht);

% Filter for positive solutions
B_ht = B_ht(B_ht > 0);
C_ht_tip = C_ht_tip(C_ht_tip > 0);
C_ht_root = C_ht_root(C_ht_root > 0);
meanAerodynamicChord_ht = meanAerodynamicChord_ht(meanAerodynamicChord_ht > 0);

% Display results
if ~isempty(B_ht) && ~isempty(C_ht_tip) && ~isempty(C_ht_root) && ~isempty(meanAerodynamicChord_ht)
    fprintf('Horizontal Tail Span (B_ht): %.2f m\n', B_ht(1));
    fprintf('Horizontal Tail Root Chord (C_ht_root): %.2f m\n', C_ht_root(1));
    fprintf('Horizontal Tail Tip Chord (C_ht_tip): %.2f m\n', C_ht_tip(1));
    fprintf('Mean Aerodynamic Chord (meanAerodynamicChord_ht): %.2f m\n', meanAerodynamicChord_ht(1));
else
    fprintf('No positive solutions found.\n');
end





% Define variables
AR_v = 1.8;                  % Aspect Ratio of the horizontal tail
taperRatio_lambdav = 0.5;    % Taper ratio of the horizontal tail
S_vt = Sv;                   % Planform area of the horizontal tail in m^2

% Solve for the mean aerodynamic chord
syms meanAerodynamicChord_vt C_vt_root B_vt C_vt_tip

% Equation for the mean aerodynamic chord
meanAerodynamicChord_eqn = (2/3) * C_vt_root * ...
    (1 + taperRatio_lambdav + taperRatio_lambdav^2) / ...
    (1 + taperRatio_lambdav) == meanAerodynamicChord_vt;

% Equation for the span of the horizontal tail
AR_eqn = AR_v == B_vt / meanAerodynamicChord_vt;

% Equation for the taper ratio
lambda_eqn = taperRatio_lambdav == C_vt_tip / C_vt_root;

% Equation for the planform area
S_vt_eqn = B_vt * meanAerodynamicChord_vt == S_vt;

% Solve the system of equations
sol = solve([meanAerodynamicChord_eqn, S_vt_eqn, AR_eqn, lambda_eqn], ...
    [B_vt, C_vt_tip, C_vt_root, meanAerodynamicChord_vt]);

% Extract solutions and filter for positive values
B_vt = double(sol.B_vt);
C_vt_tip = double(sol.C_vt_tip);
C_vt_root = double(sol.C_vt_root);
meanAerodynamicChord_vt = double(sol.meanAerodynamicChord_vt);

% Filter for positive solutions
B_vt = B_vt(B_vt > 0);
C_vt_tip = C_vt_tip(C_vt_tip > 0);
C_vt_root = C_vt_root(C_vt_root > 0);
meanAerodynamicChord_vt = meanAerodynamicChord_vt(meanAerodynamicChord_vt > 0);

% Display results
if ~isempty(B_vt) && ~isempty(C_vt_tip) && ~isempty(C_vt_root) && ~isempty(meanAerodynamicChord_vt)
    fprintf('Vertical Tail Span (B_vt): %.2f m\n', B_vt(1));
    fprintf('Vertical Tail Root Chord (C_vt_root): %.2f m\n', C_vt_root(1));
    fprintf('Vertical Tail Tip Chord (C_vt_tip): %.2f m\n', C_vt_tip(1));
    fprintf('Mean Aerodynamic Chord (meanAerodynamicChord_vt): %.2f m\n', meanAerodynamicChord_vt(1));
else
    fprintf('No positive solutions found.\n');
end

