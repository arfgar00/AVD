function optIndex = wingMDO(x)
    global airfoil cruise1 cruise2 C_cruise Wto Swet_SrefWing Swet_SrefBody l_d bdes bodyDiameter
    
    %x is the array of wing geometry to be optimized
    %x = [cr ck ct Lambdain50 Lambdaout50 yk twistmax]
    wing = WingGeometry();
    wing.cr = x(1);
    wing.ck = x(2);
    wing.ct = x(3);
    wing.s = bdes/2;
    wing.Lambdain50 = x(4);
    wing.Lambdaout50 = x(5);
    wing.yk = x(6);
    wing.twist_max = x(7);
    wing = wing.calcSref();
    %disp(["bdes", bdes])
    wing.N = 501;
    wing = wing.createStrips();
    wing = wing.calcSc();

    cruise1 = cruise1.init(wing.cbar);
    cruise2 = cruise2.init(wing.cbar);

    [CLclean1, CDi1, Cly1] = LLESwept(wing,airfoil,cruise1);
    CDF1 = CDFfun(cruise1.M, cruise1.Re, wing, airfoil.x_cm, airfoil.t_c, l_d, Swet_SrefWing, Swet_SrefBody);
    CDW1 = CDWfun(wing,Cly1,cruise1.M,airfoil.t_c);
    CDtotal1 = CDF1 + CDW1 + CDi1;
    L_D1 = CLclean1/CDtotal1;
    %disp(["L/D cruise 1 = " L_D1])

    [CLclean2, CDi2, Cly2] = LLESwept(wing,airfoil,cruise2);
    CDF2 = CDFfun(cruise2.M, cruise2.Re, wing, airfoil.x_cm, airfoil.t_c, l_d, Swet_SrefWing, Swet_SrefBody);
    CDW2 = CDWfun(wing,Cly2,cruise2.M,airfoil.t_c);
    CDtotal2 = CDF2 + CDW2 + CDi2;
    L_D2 = CLclean2/CDtotal2;
    %disp(["L/D cruise 2 = " L_D2])
    figure(2)
    clf;
    wing.plotWing()
    
    %Structural weight
    Ww = wingWeightLTH(wing,airfoil.t_c,Wto);

    %Performance
    W_f_0 = Wf_W0(cruise1, cruise2, L_D1, L_D2, C_cruise);
    optIndex = Ww*L_D1; %index to minimize
end