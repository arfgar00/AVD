function optIndex = wingMDO(x)
    %airfoil name:
    polarName = "xf-sc21010-il-1000000.csv";
    shapeName = "sc21010.dat.txt";
    
    %x is the array of wing geometry to be optimized
    %x = [cr ck ct s Lambdain50 Lambdaout50 yk]
    
    %Cruise condition is hard-coded in this function:
    cruise = AirCondition();
    cruise.M = 0.83;
    cruise.h = convlength(39000, 'ft','m');
    h_cruise2 = convlength(18000, 'ft','m');
    v_of_second_cruise = convvel(400, 'kts', 'm/s');
    C_cruise = 0.5 ./ 3600;

    %Take off weight hard coded
    Wto = 4.2188*1e5*9.81;

    %Wetted area ratio is hard-coded in this function:
    Swet_SrefWing = 2.2;
    Swet_SrefBody = 5;
    l_d = 10;

    airfoil = Airfoil();
    airfoil = airfoil.readPolar(polarName);
    airfoil = airfoil.readShape(shapeName);
    airfoil = airfoil.interpShape(9);
    
    wing = WingGeometry();
    wing.cr = 12;
    wing.ck = x(1);
    wing.ct = x(2);
    wing.s = 65/2;
    wing.Lambdain50 = x(3);
    wing.Lambdaout50 = x(4);
    wing.yk = x(5);
    wing.twist_max = x(6);
    wing = wing.calcSref();
    wing.N = 50;
    wing = wing.createStrips();
    wing = wing.calcSc();

    cruise = cruise.init(wing.cbar);

    [CLclean, CDi, Cly] = LLE_new(wing,airfoil,cruise,40,2e-1);
    CDF = CDFfun(cruise.M, cruise.Re, wing, airfoil.x_cm, airfoil.t_c, l_d, Swet_SrefWing, Swet_SrefBody);
    CDW = CDWfun(wing,Cly,cruise.M,airfoil.t_c);
    CDtotal = CDF + CDW + CDi;
    disp(["CDF = ", num2str(CDF)])
    disp(["CDW = ", num2str(CDW)])
    disp(["CDi = ", num2str(CDi)])
    disp(["CDtotal = ", num2str(CDtotal)])
    disp(["CLclean = ", num2str(CLclean)])
    L_D = CLclean/CDtotal;
    disp(['L/D = ', num2str(L_D)]);
    disp(["Sref = ",num2str(wing.SREF)])
    Ww = wingWeightLTH(wing,airfoil.t_c,Wto);
    
    [T_65, c_65, P_65, rho_65] = atmosisa(h_cruise2);
    M_cruise2 =  v_of_second_cruise/c_65;
    %[W_f_0, W_i] = Wf_W0(cruise.h, h_cruise2, cruise.M, M_cruise2, L_D, C_cruise)
    optIndex = Ww/L_D %index to minimize
end