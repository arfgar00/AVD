function T_T0 = T_T0(h)
    rho = @(h) 1.225*ISApressure(h)/ISApressure(0)*ISAtemp(0)/ISAtemp(h);
    rho0 = rho(0);
    rho1 = rho(11e3);
    if h < 11e3
        T_T0 = (rho(h)/rho0)^0.7;
    else
        T_T0 = (rho(11e3)/rho0)^(-0.3)*rho(h)/rho(0);
    end
end