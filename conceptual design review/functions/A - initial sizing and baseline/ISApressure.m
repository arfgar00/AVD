function p = ISApressure(h)
    h1 = 11e3;
    p1 = 22644;
    p0 = 101325;
    T0 = 288.2;
    g0 = 9.81;
    lambda0 = -0.0065;
    R = 287;
    T1 = ISAtemp(h1);
    if h < h1
        p = p0*(1 + lambda0*h/T0)^(-g0/(R*lambda0));
    else
        p = p1*exp(-g0/(R*T1)*(h - h1));
    end
end