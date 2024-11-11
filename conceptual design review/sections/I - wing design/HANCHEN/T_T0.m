function T_T0 = T_T0(h)
% This function is used to calculate thrust of a jet engine
% at different altitude by ignoring all the effects of 
% airspeed, which is applicable above around 30000ft with
% Mach number larger than 0.6
%
% T_T0 = T_T0(h)
%
% outputs: T_T0 - Thrust ratio
% inputs:  h    - altitude (m)

[T0, c0, P0, rho0] = atmosisa(0);
[T1, c1, P1, rho1] = atmosisa(11e3);
[Th, ch, Ph, rhoh] = atmosisa(h);

if h < 11e3
    T_T0 = (rhoh / rho0) ^ 0.7;
else
    T_T0 = (rho1 / rho0) ^ (-0.3) * (rhoh / rho0);
end

end