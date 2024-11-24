function x_ac = x_ac_M(x_ac_0, M, S_wing)
% This function is used to calculate the shift of lift surface
% aerodynamic center due to the compressible effect of flow
%
% x_ac_M(x_ac_0, M, S_w)
%
% outputs: x_ac      - Aerodynamic center at M=0 (m)
% inputs:  x_ac_0    - Aerodynamic center at M=0 (m)
%          M         - Mach number
%          S_wing    - Reference wing area (m^2)

if M < 0

    error('Mach number cannot smaller than 0')
    
elseif M < 0.4

    x_ac = x_ac_0;

elseif M <= 1.1

    x_ac = x_ac_0 + 0.26 * ((M - 0.4) ^ 2.5) * ((S_wing) ^ 0.5);

else

    error('Mach number larger than 1.1, equation invalid')

end

end