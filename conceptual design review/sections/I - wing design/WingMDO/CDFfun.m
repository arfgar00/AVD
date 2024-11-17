function CDFnum = CDFfun(M, Re, mywing, x_cm, t_c, l_d, Swet_Srefwing, Swet_Srefbody)
    %Schlichting
    CF = (2*log10(Re)- 0.65)^(-2.3);
    %Eqn.6, lambdam is the sweep angle at maximum thickness chord position
    FFWing = @(Lambdam, x_cm, t_c, M) (1 + 0.6./x_cm .* t_c + 100*t_c.^4) *(1.34*M.^0.18*cos(Lambdam).^(0.28));
    %Eqn.10
    FFBody = @(l_d) 1 + 0.0025 * l_d + 60 / (l_d)^3;
    %Use weighted avraged FFWing between outer and inner sections
    FFWingval = mywing.Sin/mywing.SREF*FFWing(mywing.Lambdax_c(x_cm,0), x_cm, t_c, M) + mywing.Sout/mywing.SREF*FFWing(mywing.Lambdax_c(x_cm,mywing.s), x_cm, t_c, M);
    FFBodyval = FFBody(l_d);
    CDFnum = CF*(FFWingval*Swet_Srefwing + FFBodyval*Swet_Srefbody);
end