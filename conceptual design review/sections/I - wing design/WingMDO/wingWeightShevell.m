function Ww = wingWeightShevell(mywing, myairfoil, nult, Wto, Wf)
    % wzf is zero fuel weight = Wto - Wf
    % Iw in imperial units
    Wzf = Wto - Wf;
    Wto = convforce(Wto, "N","lbf");
    Wzf = convforce(Wzf, "N","lbf");
   
    lambda = mywing.ct/mywing.cr;
    b = convlength(mywing.b, "m", "ft");
    t_c = convlength(myairfoil.t_c, "m", "ft");
    Sref = mywing.SREF*10.7639;
    Iw = nult*b^3*sqrt(Wto * Wzf)/(t_c * Sref * cos(mywing.Lambdaout50)^2) * (1+2*lambda)/(1 + lambda) * 1e-6;
    Ww = mywing.SREF* (4.22 + 1.642 * Iw);
end