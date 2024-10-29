function T = ISAtemp(h)
    if h < 11e3
        T = 288.15 + -0.0065*h;  % Temperature decreases with altitude below 11 km
    else
        T = 216.65;  % Temperature is constant at 216.65 K above 11 km
    end
end
