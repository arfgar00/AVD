function optIndex = wingMDO(x)
    %global airfoil cruise1 cruise2 C_cruise Wto Swet_SrefWing Swet_SrefBody l_d bdes W_ZF_Nowing Ww nult Wpayload bodyLength bodyDiameter
    run("SetParametersMDO.m") %ensure functions are independent

    %x is the array of wing geometry to be optimized
    %x = [cr ck ct Lambdain50 Lambdaout50 yk twistmax]
    wing = x2wing(x);
    % wing = WingGeometry();
    % wing.cr = x(1);
    % wing.ck = x(2);
    % wing.ct = x(3);
    % wing.s = bdes/2;
    % wing.Lambdain50 = x(4);
    % wing.Lambdaout50 = x(5);
    % wing.yk = x(6);
    % wing.twist_max = x(7);
    % wing = wing.calcSref();
    % %disp(["bdes", bdes])
    % wing.N = 301;
    % wing = wing.createStrips();
    % wing = wing.calcSc();

    cruise1 = cruise1.init(wing.cbar);
    cruise2 = cruise2.init(wing.cbar);
        
    Swet_SrefBody = pi*bodyDiameter*bodyLength*1.1/wing.SREF;
    
    [CLclean1, CDi1, Cly1, CDprofile1] = LLESwept(wing,airfoil,cruise1, Damping, tolerance, bodyDiameter);
    CDF1 = CDFfun(cruise1.M, cruise1.Re, wing, airfoil.x_cm, airfoil.t_c, l_d, Swet_SrefWing, Swet_SrefBody);
    CDW1 = CDWfun(wing,Cly1,cruise1.M,airfoil.t_c);
    CDtotal1 = CDF1 + CDW1 + CDi1 + CDprofile1;
    L_D1 = CLclean1/CDtotal1;
    disp(["L/D cruise 1 = " L_D1])

    [CLclean2, CDi2, Cly2, CDprofile2] = LLESwept(wing,airfoil,cruise2, Damping, tolerance, bodyDiameter);
    CDF2 = CDFfun(cruise2.M, cruise2.Re, wing, airfoil.x_cm, airfoil.t_c, l_d, Swet_SrefWing, Swet_SrefBody);
    CDW2 = CDWfun(wing,Cly2,cruise2.M,airfoil.t_c);
    CDtotal2 = CDF2 + CDW2 + CDi2 + CDprofile2;
    L_D2 = CLclean2/CDtotal2;
    disp(["L/D cruise 2 = " L_D2])
    %figure(2)
    %clf;
    %wing.plotWing()
    
    %Performance
    W_f_0 = Wf_W0(cruise1, cruise2, L_D1, L_D2, C_cruise);
    disp(["Fuel fraction: ", W_f_0])
    %Structural weight
    Ww = wingWeightRaymer(wing,airfoil,Wto0,nult,wing.SREF*0.0612);

    %update total weight
    Wto = (W_ZF_Nowing + Ww + Wpayload)/(1 - W_f_0);
    %disp(["L/D1", ""])
    optIndex = Wto; %index to minimize
    disp(["Wing Mass [kg]= ", Ww/9.81])
    disp(["Total Mass [kg] = ", Wto/9.81])
end